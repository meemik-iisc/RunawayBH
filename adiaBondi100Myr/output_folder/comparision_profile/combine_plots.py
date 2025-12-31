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

FOLDERS = ["dens", "pres", "mass_flux", "temp", "velr", "time_scales"]
Plot_name = ["Density", "Pressure","Mass Flux", "Temperature", "Radial Velocity", "Time Scales"]
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
import os
import cv2
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image
from pathlib import Path

def create_subplots():
    """Create 2x3 high-resolution subplots from 6 folders and save to output folder."""
    
    # Create output folder if it doesn't exist
    output_path = Path(OUTPUT_FOLDER)
    output_path.mkdir(exist_ok=True)
    print(f"✓ Created output folder: {OUTPUT_FOLDER}")
    
    # Get reference image size from first folder, first file
    sample_path = Path(FOLDERS[0]) / "t_0_comparison.png"
    if not sample_path.exists():
        print(f"✗ Error: Sample image {sample_path} not found!")
        return False
    
    sample_img = Image.open(sample_path)
    img_w, img_h = sample_img.size  # Original PNG pixel dimensions
    
    # Layout: 3 columns, 2 rows
    cols, rows = 3, 2
    target_w = cols * img_w
    target_h = rows * img_h
    
    # High DPI for sharp output
    dpi = 300
    fig_w_inch = target_w / dpi
    fig_h_inch = target_h / dpi
    
    print(f"✓ Input image size: {img_w}x{img_h}")
    print(f"✓ Target canvas: {target_w}x{target_h} pixels at {dpi} DPI")
    
    # Loop through each time step
    for t in range(NUM_FILES):
        output_file_name = f"t_{t}.png"
        
        # Create high-res figure matching input pixel budget
        fig, axes = plt.subplots(rows, cols, figsize=(fig_w_inch, fig_h_inch), dpi=dpi)
        axes = axes.flatten()
        
        # Extract time number for super title
        time_value = t * time_scale
        fig.suptitle(f"Time: {time_value} {time_label}", fontsize=16, y=0.98)
        
        # Load and display image from each folder
        for idx, folder in enumerate(FOLDERS):
            filename = f"t_{t}_comparison.png"
            image_path = Path(folder) / filename
            
            try:
                # Check if file exists
                if not image_path.exists():
                    print(f"⚠ Warning: {image_path} not found")
                    axes[idx].text(0.5, 0.5, f"File not found\n{filename}", 
                                  ha='center', va='center', fontsize=12, transform=axes[idx].transAxes)
                    axes[idx].set_title(Plot_name[idx], fontsize=14, pad=10)
                    axes[idx].axis('off')
                    continue
                
                # Load image using PIL (preserves original resolution)
                image = Image.open(image_path)
                
                # Display WITHOUT interpolation for pixel-perfect sharpness
                axes[idx].imshow(image, interpolation='nearest')
                axes[idx].set_title(Plot_name[idx], fontsize=14, pad=10)
                axes[idx].axis('off')  # Hide axes
                
            except Exception as e:
                print(f"✗ Error loading {image_path}: {e}")
                axes[idx].text(0.5, 0.5, f"Error loading\n{filename}", 
                              ha='center', va='center', fontsize=12, transform=axes[idx].transAxes)
                axes[idx].set_title(Plot_name[idx], fontsize=14, pad=10)
                axes[idx].axis('off')
        
        # Minimal tight layout to avoid overlap but preserve space
        plt.tight_layout(pad=0.5)
        
        # Save at high DPI with full quality
        output_file = output_path / output_file_name
        plt.savefig(output_file, dpi=dpi, bbox_inches='tight', 
                   facecolor='white', edgecolor='none')
        print(f"✓ Saved: {output_file} ({target_w}x{target_h} px)")
        
        # Close figure to free memory
        plt.close(fig)
    
    print(f"\n✓ All {NUM_FILES} high-res files created successfully!")
    return True

def create_video_from_subplots():
    """Create high-quality MP4 video from PNG files using OpenCV."""
    
    print("=" * 70)
    print("STEP 2: CREATING HIGH-QUALITY MP4 VIDEO")
    print("=" * 70)
    print()
    
    data_path = Path(OUTPUT_FOLDER)
    
    if not data_path.exists():
        print(f"✗ Error: Folder '{OUTPUT_FOLDER}' not found!")
        return False
    
    # Get all PNG files sorted NUMERICALLY
    png_files = sorted(data_path.glob("t_*.png"), key=lambda x: int(x.stem.split('_')[1]))
    
    if not png_files:
        print(f"✗ Error: No PNG files found in '{OUTPUT_FOLDER}'!")
        return False
    
    print(f"✓ Found {len(png_files)} PNG files")
    print(f"  Files: {png_files[0].name} to {png_files[-1].name}")
    
    # Read first image to get dimensions (BGR format)
    first_image = cv2.imread(str(png_files[0]))
    if first_image is None:
        print(f"✗ Error: Could not read first image {png_files[0]}")
        return False
    
    height, width = first_image.shape[:2]
    print(f"✓ Frame resolution: {width}x{height}")
    
    # Better codec options for quality (H.264 > mp4v)
    fourcc = cv2.VideoWriter_fourcc(*'H264')  # or 'mp4v' as fallback
    video_path = Path(OUTPUT_VIDEO)
    video_path.parent.mkdir(exist_ok=True, parents=True)
    
    out = cv2.VideoWriter(str(video_path), fourcc, FPS, (width, height))
    
    if not out.isOpened():
        print("✗ H.264 failed, trying mp4v...")
        fourcc = cv2.VideoWriter_fourcc(*CODEC)  # Fallback to mp4v
        out = cv2.VideoWriter(str(video_path), fourcc, FPS, (width, height))
        if not out.isOpened():
            print(f"✗ Error: Could not open video writer with {CODEC}!")
            print("  Install FFmpeg: sudo apt install ffmpeg (Linux)")
            return False
    
    print(f"✓ Codec: {chr(fourcc&0xFF)}{chr((fourcc>>8)&0xFF)}{chr((fourcc>>16)&0xFF)}{chr((fourcc>>24)&0xFF)}")
    print(f"✓ FPS: {FPS} | Duration: ~{len(png_files)/FPS:.1f}s")
    print()
    
    # Write frames
    for idx, image_path in enumerate(png_files):
        frame = cv2.imread(str(image_path))
        if frame is None:
            print(f"⚠ Skipping {image_path.name}")
            continue
        
        out.write(frame)
        if (idx + 1) % 5 == 0 or idx == len(png_files) - 1:
            print(f"✓ Progress: {idx + 1}/{len(png_files)} ({100*(idx+1)/len(png_files):.0f}%)")
    
    out.release()
    cv2.destroyAllWindows()
    
    print(f"\n✓ Video saved: {OUTPUT_VIDEO} ({width}x{height}, {len(png_files)} frames)")
    return True
def convert_images_to_pdf():
    """Convert all PNG images from data folder to PDF with natural sorting."""
    
    print("=" * 70)
    print("STEP 3: CONVERTING IMAGES TO PDF")
    print("=" * 70)
    print()
    
    data_path = Path(OUTPUT_FOLDER)
    
    # Check if data folder exists
    if not data_path.exists():
        print(f"✗ Error: Folder '{OUTPUT_FOLDER}' not found!")
        return False
    
    # Get all PNG files with NATURAL SORT (t_0, t_1, ..., t_10, not t_0, t_1, t_10, t_2...)
    try:
        from natsort import natsorted
        png_files = natsorted(data_path.glob("t_*.png"))
    except ImportError:
        print("⚠ natsort not installed. Using basic numeric sort...")
        png_files = sorted(data_path.glob("t_*.png"), key=lambda x: int(x.stem.split('_')[1]))
    
    if not png_files:
        print(f"✗ Error: No PNG files found in '{OUTPUT_FOLDER}'!")
        return False
    
    print(f"✓ Found {len(png_files)} PNG files")
    print(f"  Files: {png_files[0].name} to {png_files[-1].name}")
    print()
    
    # Convert images to PDF
    try:
        import img2pdf
    except ImportError:
        print("✗ Error: img2pdf not installed!")
        print("  Install with: pip install img2pdf")
        return False
    
    # Create list of image paths as strings
    image_list = [str(img) for img in png_files]
    
    # Convert to PDF
    try:
        pdf_path = f"{OUTPUT_FOLDER}.pdf"
        pdf_bytes = img2pdf.convert(image_list)
        
        with open(pdf_path, "wb") as f:
            f.write(pdf_bytes)
        
        print(f"✓ PDF created successfully: {pdf_path}")
        print(f"✓ Total pages: {len(png_files)}")
        print()
        return True
        
    except Exception as e:
        print(f"✗ Error converting to PDF: {e}")
        return False

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
    
    # Step 3: Create PDF
    success_pdf = convert_images_to_pdf()
    
    if success_subplots and (success_video or success_pdf):
        print("╔" + "=" * 68 + "╗")
        print("║" + " " * 68 + "║")
        print("║" + "✓ SUCCESS: ALL TASKS COMPLETED!".center(68) + "║")
        print("║" + " " * 68 + "║")
        print("╚" + "=" * 68 + "╝")
        print()
        print(f"Output files:")
        print(f"  • Subplots folder: {OUTPUT_FOLDER}/")
        print(f"  • Video file: {OUTPUT_VIDEO}")
        print(f"  • PDF file: {OUTPUT_FOLDER}.pdf")
        print()
    else:
        print("⚠ Some tasks failed. Check output above.")
        print()
