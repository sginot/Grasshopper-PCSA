// End result wanted
// Make 2d3d convex hull and use it to mask the original greyscale images
// Greyscale images are then binarized using the different methods available in Fiji
// 3D volume is computed based on area of muscle*voxel depth
// Results stored in a table with individual specimens as lines, and binarization method as columns

input_seg = getDirectory("Input Directory for segmented muscles");
	// input directory containing all individual segmented muscle volumes (in our case MHD format)
	
input_orig = getDirectory("Input Directory for original grayscale images");
	// input directory containing original greyscale images
	
output_2d3d = getDirectory("Output directory to store 2d3d convex hulls");
	
output_results = getDirectory("Output Directory for tables containing results");

a1 = newArray("Default", "Huang", "Intermodes", "IJ_IsoData", "MaxEntropy", "Minimum", "Moments", "Otsu", "Yen");

Dialog.create("Binarization method");
Dialog.addChoice("Algorithm", a1);
Dialog.show();
binar = Dialog.getChoice();

lsfil = getFileList(input_seg);

for (i = 0; i < lsfil.length; i++) {

	if (endsWith(lsfil[i], "mhd")) {
		print("MHD Image");
		print(lsfil[i]);
		run("MHD/MHA...", "open=["+input_seg +lsfil[i]+"]");
	} else if (endsWith(lsfil[i], "raw")) {
		print("Raws should correspond to an MHD/MHA file, to be automatically opened");
		continue;
	} else {
		print("Not an appropriate image format");
		continue;
	}
	
	Dialog.create(lsfil[i]);
	Dialog.addNumber("Voxel size (mm)", "");
	Dialog.show();
	width = Dialog.getNumber();	
	run("Set Scale...", "distance=1 known="+width+" pixel=1 unit=mm");	
		
	for (j = 1; j <= nSlices; j++) {
    	setSlice(j);
    	
    	setThreshold(129, 255, "raw");
    	run("Set Measurements...", "mean redirect=None decimal=3");
    	run("Select All");
    	run("Measure");
    	mean = getResult("Mean");
    	if (mean == "0") {
    		continue;
    	} 
    	else {
   		run("Create Selection");
   		run("Convex Hull");
   		run("Clear", "slice");
  		}
	}
	
	nam = lsfil[i];
	run("MHD/MHA ...", "save=["+output_2d3d+ nam+ ".mha]");
	close("*");
	close("Results");
}

lsfil = getFileList(output_2d3d);

for (i = 0; i < lsfil.length; i++) {

	run("MHD/MHA...", "open=["+output_2d3d +lsfil[i]+"]");
	rename("mask");
	
	print("Select gray value image corresponding to "+lsfil[i]+" to be masked");
	open(""); // Select file of specimen corresponding to the opened mask
	rename("grey");
		
	imageCalculator("AND create stack", "grey","mask");
	rename("masked");
	
	Dialog.create(lsfil[i]);
	Dialog.addNumber("Voxel size (mm)", "");
	Dialog.show();
	width = Dialog.getNumber();	
	
	//binar = "Default"
	//width = 0.01198
	run("Make Binary", "method="+binar+" background=Dark calculate");
	run("Analyze Particles...", "size=20-Infinity pixel show=Masks stack");
	run("Set Scale...", "distance=1 known="+width+" pixel=1 unit=mm");
	run("Set Measurements...", "mean area redirect=None decimal=3");
	
		for (j = 1; j <= nSlices; j++) {
    	setSlice(j);
    	
    	setThreshold(129, 255, "raw");
   		run("Create Selection");
     	run("Measure");
  		}
  		
	nam = lsfil[i];
	saveAs("Results", output_results+ nam+ ".csv");	
	close("*");
	close("Results");	
}
		

	