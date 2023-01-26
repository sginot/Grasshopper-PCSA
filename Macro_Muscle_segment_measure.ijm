// End results wanted:
// Define an input folder containing raws, image sequences or else
// Pop up menu, asking for folder and file type?
// Define scale(s) and dimensions, also in popup menu
// Define  Output folder
// Run macro automatically through scans remove cuticle 
// Save that to output folder
// Run second step macro 2dto3d convex hull /!\ different macro!!!

arr = newArray("No", "Yes");
arrform = newArray("Raw Data", "Image Sequence", "MHD/MHA");

Dialog.create("Start");
Dialog.addChoice("All images have the same voxel dimensions", arr);
Dialog.show();
dimsame = Dialog.getChoice();
print(dimsame);

if (dimsame == "No") {
	print("Warning: Voxel dimensions will have to be set individually");
}
else {
	input = getDirectory("Input Directory");
	output = getDirectory("Output Directory");
	Dialog.create("Output format");
	Dialog.addChoice("Choose a format", arrform);
	Dialog.show();
	outformat = Dialog.getChoice();

	if (outformat == "Raw Data") {
		ext = ".raw";
	}
	if (outformat == "Image Sequence") {
		ext = ".tif";
	}
	if (outformat == "MHD/MHA") {
		ext = ".mhd";
	}
	
	Dialog.create("Voxel dimensions");
	Dialog.addNumber("Voxel Width", "");
	Dialog.addNumber("Voxel Height", "");
	Dialog.addNumber("Voxel Depth", "");
	Dialog.show();
	width = Dialog.getNumber();
	height = Dialog.getNumber();
	depth = Dialog.getNumber();
	print("Voxel dimensions:", width, height, depth);

	if (width == height) {
		pixratio = 1;
	}
	else {
		pixratio = width/height;
	}

	print("Pixel ratio:" +pixratio);

	lsfil = getFileList(input);
	Array.show(lsfil);
	print(lsfil.length);

	for (i = 0; i < lsfil.length; i++) {

		if (endsWith(lsfil[i], "mhd")) {
			print("MHD Image");
			run("MHD/MHA...", "open=["+input +lsfil[i]+"]");
		} else if (endsWith(lsfil[i], "mha")) {
			print("MHA Image");
			run("MHD/MHA...", "open=["+input +lsfil[i]+"]");
		} else if (endsWith(lsfil[i], "tif")) {
			print("TIF Stack");
			run("TIFF Virtual Stack...", "open=["+input +lsfil[i]+"]");
			run("Duplicate...", "duplicate");
		} else if (endsWith(lsfil[i], "tiff")) {
			print("TIFF Stack");
			run("TIFF Virtual Stack...", "open=["+input +lsfil[i]+"]");
			rename("Virtual")
			run("Duplicate...", "duplicate");
			close("Virtual");
		} else if (endsWith(lsfil[i], "/")) {
			print("Directory containing Image Stack");
			newdir = input +lsfil[i];
			newlist= getFileList(newdir);
			first = newlist[1];
			run("Image Sequence...", "open=["+newdir +first+"] sort");
		} else if (endsWith(lsfil[i], "raw")) {
			print("Raws should correspond to an MHD/MHA file, to be automatically opened");
			continue;
		} else {
			print("Not an appropriate image format");
		}

// Actual image treatment starts here		
		run("Set Scale...", "distance=1 known="+width+" pixel="+pixratio+" unit=mm");
		run("Set Measurements...", "area redirect=None decimal=3");
		run("Make Binary", "method=Default background=Default calculate");
		rename("original");
		// Carefully choose the size value in the next line
		run("Analyze Particles...", "size=0.20-Infinity show=Masks stack");
		rename("mask");
		imageCalculator("XOR create stack", "original","mask");
		rename("masked");
		selectWindow("masked");

		if (endsWith(lsfil[i], "/")) {
			nam = lsfil[i];
			newnam = substring(nam, 0, nam.length-1);
		}
		else {
			nam = lsfil[i];
			newnam = substring(nam, 0, nam.length-4);
		}

		if (outformat == "Image Sequence") {
			out = output +newnam+ "/";
			print(out);
			print(out +newnam+ "0000" +ext);
			rename(newnam);
			File.makeDirectory(out);
			run("Image Sequence... ", "format=TIFF save=[" +out +newnam+ "0000" +ext+ "]");
		}
		if (outformat == "Raw Data") {
			print(output+ newnam+ ext);
			saveAs("Raw Data", output+ newnam+ ext);
		}
		if (outformat == "MHD/MHA") {
			print(output+ newnam+ ext);
			run("MHD/MHA ...", "save=[" +output+ newnam+ ext+ "]");
		}
		
		close("*");
	}
}


