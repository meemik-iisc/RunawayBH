"""
Script to create subplots from multiple folders and save as PNG files.

This script reads PNG files from 4 folders (dens, pres, temp, velr),
creates 2x2 subplots for each time step, and saves them to a 'data' folder.

Folder structure:
    dens/
        ├── t_0.png
        ├── t_1.png
        └── ...
    pres/
    temp/
    velr/
    data/  (output folder - will be created)
        ├── t_0.png
        ├── t_1.png
        └── ...
"""

import os
import cv2
import matplotlib.pyplot as plt
from PIL import Image
from pathlib import Path

# ============================================================================
# CONFIGURATION
# ============================================================================

FOLDERS = ["dens", "pres", "temp", "velr"]
Plot_name = ["Density", "Pressure", "Temperature", "Radial Velocity"]
OUTPUT_FOLDER = "data"
OUTPUT_VIDEO = "data_video.mp4"
NUM_FILES = 11
FPS = 2  # Frames per second for video
CODEC = "mp4v"  # Video codec
time_scale=10.0
time_label = "Myr"
# ============================================================================
# MAIN SCRIPT
# ============================================================================

def create_subplots():
    """Create 2x2 subplots from 4 folders and save to output folder."""
    
    # Create output folder if it doesn't exist
    output_path = Path(OUTPUT_FOLDER)
    output_path.mkdir(exist_ok=True)
    print(f"✓ Created output folder: {OUTPUT_FOLDER}")
    
    # Loop through each time step
    for t in range(NUM_FILES):
        filename = f"t_{t}.png"
        
        # Create figure with 2x2 subplots
        fig, axes = plt.subplots(2, 2, figsize=(12, 10))
        axes = axes.flatten()
        
        # Extract time number from filename for super title
        time_value = t  # You can modify this if needed
        fig.suptitle(f"Time: {time_value*time_scale} {time_label}")
        
        # Load and display image from each folder
        for idx, folder in enumerate(FOLDERS):
            try:
                # Construct image path
                image_path = Path(folder) / filename
                
                # Check if file exists
                if not image_path.exists():
                    print(f"⚠ Warning: {image_path} not found")
                    axes[idx].text(0.5, 0.5, f"File not found\n{filename}", 
                                  ha='center', va='center', fontsize=12)
                    axes[idx].set_title(Plot_name[idx])
                    axes[idx].axis('off')
                    continue
                
                # Load image using PIL
                image = Image.open(image_path)
                
                # Display image in subplot
                axes[idx].imshow(image)
                axes[idx].set_title(Plot_name[idx])
                axes[idx].axis('off')  # Hide axes
                
            except Exception as e:
                print(f"✗ Error loading {image_path}: {e}")
                axes[idx].text(0.5, 0.5, f"Error loading\n{filename}", 
                              ha='center', va='center', fontsize=12)
                axes[idx].set_title(Plot_name[idx])
                axes[idx].axis('off')
        
        # Adjust layout to prevent overlap
        plt.tight_layout()
        
        # Save figure
        output_file = output_path / filename
        plt.savefig(output_file, dpi=100, bbox_inches='tight')
        print(f"✓ Saved: {output_file}")
        
        # Close figure to free memory
        plt.close(fig)
    
    print(f"\n✓ All {NUM_FILES} files created successfully in '{OUTPUT_FOLDER}' folder!")
    print()
    return True
    
def create_video_from_subplots():
    """Create MP4 video from PNG files using OpenCV."""
    
    print("=" * 70)
    print("STEP 2: CREATING MP4 VIDEO")
    print("=" * 70)
    print()
    
    data_path = Path(OUTPUT_FOLDER)
    
    # Check if data folder exists
    if not data_path.exists():
        print(f"✗ Error: Folder '{OUTPUT_FOLDER}' not found!")
        return False
    
    # Get all PNG files sorted NUMERICALLY by the number in filename
    png_files = sorted(data_path.glob("t_*.png"), key=lambda x: int(x.stem.split('_')[1]))
    
    if not png_files:
        print(f"✗ Error: No PNG files found in '{OUTPUT_FOLDER}' folder!")
        return False
    
    print(f"✓ Found {len(png_files)} PNG files")
    print(f"  Files: {png_files[0].name} to {png_files[-1].name}")
    
    # Read first image to get dimensions
    first_image = cv2.imread(str(png_files[0]))
    if first_image is None:
        print(f"✗ Error: Could not read first image {png_files[0]}")
        return False
    
    height, width = first_image.shape[:2]
    print(f"✓ Image dimensions: {width}x{height}")
    
    # Create video writer
    fourcc = cv2.VideoWriter_fourcc(*CODEC)
    out = cv2.VideoWriter(OUTPUT_VIDEO, fourcc, FPS, (width, height))
    
    if not out.isOpened():
        print(f"✗ Error: Could not open video writer!")
        print("  Make sure FFmpeg is installed on your system")
        return False
    
    print(f"✓ Video codec: {CODEC}")
    print(f"✓ Video FPS: {FPS}")
    print(f"\nCreating video...\n")
    
    # Write images to video
    for idx, image_path in enumerate(png_files):
        try:
            frame = cv2.imread(str(image_path))
            if frame is None:
                print(f"⚠ Warning: Could not read {image_path.name}")
                continue
            
            out.write(frame)
            print(f"✓ Frame {idx + 1}/{len(png_files)}: {image_path.name}")
            
        except Exception as e:
            print(f"✗ Error processing {image_path.name}: {e}")
    
    # Release video writer
    out.release()
    
    print(f"\n✓ Video created successfully: {OUTPUT_VIDEO}")
    print(f"✓ Total frames: {len(png_files)}")
    print(f"✓ Video duration: ~{len(png_files) / FPS:.2f} seconds")
    print()
    
    return True



if __name__ == "__main__":
    print("\n")
    print("╔" + "=" * 68 + "╗")
    print("║" + " " * 68 + "║")
    print("║" + "  SUBPLOT GENERATOR + VIDEO CREATOR".center(68) + "║")
    print("║" + " " * 68 + "║")
    print("╚" + "=" * 68 + "╝")
    print()
    
    # Verify all source folders exist
    print("Verifying source folders...")
    all_folders_exist = True
    for folder in FOLDERS:
        if not Path(folder).exists():
            print(f"✗ Error: Folder '{folder}' not found!")
            all_folders_exist = False
        else:
            print(f"✓ Found folder: {folder}")
    
    if not all_folders_exist:
        print("\n✗ FAILED: Not all source folders found!")
        exit(1)
    
    print()
    
    # Step 1: Create subplots
    success_subplots = create_subplots()
    
    if not success_subplots:
        print("✗ FAILED: Could not create subplots!")
        exit(1)
    
    # Step 2: Create video
    success_video = create_video_from_subplots()
    
    if success_video:
        print("╔" + "=" * 68 + "╗")
        print("║" + " " * 68 + "║")
        print("║" + "✓ SUCCESS: All tasks completed!".center(68) + "║")
        print("║" + " " * 68 + "║")
        print("╚" + "=" * 68 + "╝")
        print()
        print(f"Output files:")
        print(f"  • Subplots folder: {OUTPUT_FOLDER}/")
        print(f"  • Video file: {OUTPUT_VIDEO}")
        print()
    else:
        print("╔" + "=" * 68 + "╗")
        print("║" + " " * 68 + "║")
        print("║" + "✗ FAILED: Could not create video".center(68) + "║")
        print("║" + " " * 68 + "║")
        print("╚" + "=" * 68 + "╝")
        print()
        print("Troubleshooting:")
        print("1. Make sure OpenCV is installed:")
        print("   pip install opencv-python")
        print("2. Make sure FFmpeg is installed:")
        print("   Ubuntu/Debian: sudo apt-get install ffmpeg")
        print("   macOS: brew install ffmpeg")
        print("   Windows: https://ffmpeg.org/download.html")
        print()
