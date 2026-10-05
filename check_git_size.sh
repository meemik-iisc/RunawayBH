#!/bin/bash
# ============================================================
# check_git_size.sh
# Checks sizes of files staged for git add . and warns if
# any file or total size exceeds GitHub's limits.
# ============================================================

# GitHub limits
MAX_FILE_MB=100        # GitHub hard limit per file
WARN_FILE_MB=50        # Warn if file exceeds this
MAX_REPO_MB=2048       # GitHub repo hard limit (2GB)
WARN_REPO_MB=1000      # Warn if total exceeds this (1GB)

RED='\033[0;31m'
YELLOW='\033[1;33m'
GREEN='\033[0;32m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

echo ""
echo "======================================"
echo "   Git Push Size Checker"
echo "======================================"
echo ""

# Check if inside a git repo
if ! git rev-parse --git-dir > /dev/null 2>&1; then
  echo -e "${RED}ERROR: Not inside a git repository!${NC}"
  exit 1
fi

# ---- 1. Check untracked + modified files (what 'git add .' would stage) ----
echo -e "${BLUE}[1] Scanning files that 'git add .' would stage...${NC}"
echo ""

PROBLEM_FILES=()
TOTAL_NEW_MB=0
FILE_COUNT=0

# Get all untracked and modified files
while IFS= read -r filepath; do
  [ -f "$filepath" ] || continue
  SIZE_BYTES=$(stat -c%s "$filepath" 2>/dev/null || stat -f%z "$filepath" 2>/dev/null)
  SIZE_MB=$(echo "scale=2; $SIZE_BYTES/1024/1024" | bc)
  FILE_COUNT=$((FILE_COUNT + 1))
  TOTAL_NEW_MB=$(echo "$TOTAL_NEW_MB + $SIZE_MB" | bc)

  if (( $(echo "$SIZE_MB > $MAX_FILE_MB" | bc -l) )); then
    echo -e "  ${RED}✗ BLOCKED  ${SIZE_MB}MB  →  $filepath${NC}"
    PROBLEM_FILES+=("BLOCKED|$SIZE_MB|$filepath")
  elif (( $(echo "$SIZE_MB > $WARN_FILE_MB" | bc -l) )); then
    echo -e "  ${YELLOW}⚠ WARNING  ${SIZE_MB}MB  →  $filepath${NC}"
    PROBLEM_FILES+=("WARNING|$SIZE_MB|$filepath")
  fi

done < <(git ls-files --others --modified --exclude-standard)

if [ $FILE_COUNT -eq 0 ]; then
  echo -e "  ${GREEN}No new/modified files to stage.${NC}"
else
  echo ""
  echo -e "  Total new files to be staged: ${FILE_COUNT}"
  echo -e "  Total size of new files: ${TOTAL_NEW_MB} MB"
fi

# ---- 2. Check already committed history size ----
echo ""
echo -e "${BLUE}[2] Checking current git history size...${NC}"
echo ""

HISTORY_MB=$(du -sm .git | awk '{print $1}')
echo -e "  Current .git size: ${HISTORY_MB} MB"

# ---- 3. Estimate total after push ----
echo ""
echo -e "${BLUE}[3] Push size estimate...${NC}"
echo ""

ESTIMATED_TOTAL=$(echo "$HISTORY_MB + $TOTAL_NEW_MB" | bc)
echo -e "  Estimated total after push: ${ESTIMATED_TOTAL} MB"

# ---- 4. Final verdict ----
echo ""
echo "======================================"
echo "   VERDICT"
echo "======================================"

BLOCKED=0
for entry in "${PROBLEM_FILES[@]}"; do
  TYPE=$(echo $entry | cut -d'|' -f1)
  [ "$TYPE" == "BLOCKED" ] && BLOCKED=1
done

if [ $BLOCKED -eq 1 ]; then
  echo -e "${RED}✗ CANNOT PUSH — files exceed GitHub's 100MB per-file limit!${NC}"
  echo ""
  echo -e "  Fix options:"
  echo -e "  • Add to .gitignore and don't track these files"
  echo -e "  • Use Git LFS:  git lfs track '*.png'"
  echo -e "  • Remove from history:  git filter-repo --path <file> --invert-paths"
elif (( $(echo "$ESTIMATED_TOTAL > $MAX_REPO_MB" | bc -l) )); then
  echo -e "${RED}✗ CANNOT PUSH — estimated repo size exceeds GitHub's 2GB limit!${NC}"
elif (( $(echo "$ESTIMATED_TOTAL > $WARN_REPO_MB" | bc -l) )); then
  echo -e "${YELLOW}⚠ CAN PUSH but repo is getting large (${ESTIMATED_TOTAL}MB). Consider cleaning up.${NC}"
elif [ ${#PROBLEM_FILES[@]} -gt 0 ]; then
  echo -e "${YELLOW}⚠ CAN PUSH but some large files detected. Review warnings above.${NC}"
else
  echo -e "${GREEN}✓ SAFE TO PUSH! Estimated size: ${ESTIMATED_TOTAL} MB${NC}"
fi

echo ""
echo "  Run 'git push origin master' when ready."
echo ""