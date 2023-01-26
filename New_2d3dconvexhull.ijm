// Analyze particles and multiply by voxel depth then sum
// Save resulting image sequence to second folder
// Save measurements to the same folder


input = getDirectory("Input Directory");
output = getDirectory("Output Directory");

lsfil = getFileList(input);

for (i = 0; i < lsfil.length; i++) {

	if (endsWith(lsfil[i], "mhd")) {
		print("MHD Image");
		run("MHD/MHA...", "open=["+input +lsfil[i]+"]");
	} else if (endsWith(lsfil[i], "raw")) {
		print("Raws should correspond to an MHD/MHA file, to be automatically opened");
		continue;
	} else {
		print("Not an appropriate image format");
		continue;
	}
	
	Dialog.create("Voxel size (mm)");
	Dialog.addNumber("Size", "");
	Dialog.show();
	width = Dialog.getNumber();
	run("Set Scale...", "distance=1 known="+width+" pixel=1 unit=mm");
	run("Invert", "stack");

	for (j = 1; j <= nSlices; j++) {
    	setSlice(j);
    
    	run("Set Measurements...", "mean redirect=None decimal=3");
    	run("Select All");
    	run("Measure");
    	mean = getResult("Mean");
    	if (mean == "255") {
    		continue;
    	} 
    	else {
   		run("Create Selection");
   		run("Convex Hull");
   		run("Clear", "slice");
  		}    
	}
	close("Results");
	
	run("Set Measurements...", "area redirect=None decimal=3");
	run("Analyze Particles...", "display stack");
	
	nam = lsfil[i];
	saveAs("Results", output+ nam+ ".csv");
	close("*");
}
