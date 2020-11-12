// Author: Socheata Ly, University of Massachusetts Medical School
// soki.ly@gmail.com
// INPUT: Folder containing multichannel images
// OUTPUT: Text file of nuclear and cytoplasmic counts and subfolders for each image analyzed


macro_name = "RNAscope_nuc_cyto_counter_v3b.ijm";
dbug = false;
save_files = true;
save_stats = true;


run("Conversions...", "scale");
// 123 for my images, 231 for Marie's images
rearrange = "231"; // enter number as a string ("123") = no change; macro expects C1=Green, C2=Red, C3=DAPI
C1_s1 = 2; // Difference of Gaussian settings (typical setings = 2,4,0,Otsu)
C1_s2 = 4;
C1_e = 0;
C2_s1 = C1_s1;
C2_s2 = C1_s2;
C2_e = C1_e;
auto_thr = "Otsu";
min_thr = 250; // default 250
max_thr = 2000; // default 750 (SL) or 2000 (MD)
min_size = 25; // for 3D object counter
DAPI_gb = 10;
DAPI_dilate = 0;
DAPI_thr = "Default";
//output_suffix = "_s"+C1_s1+"_s"+C1_s2+"_e"+C1_e+"_minSize"+min_size+"_thr"+auto_thr+"-"+min_thr+"-"+max_thr;
output_suffix = "_DOG"+C1_s1+"-"+C1_s2+"_DAPIdilate"+DAPI_dilate+"_minSize"+min_size+"_thr"+auto_thr+"-"+min_thr+"-"+max_thr;
output = "RNAscope_nuc_cyto_counter"+output_suffix;
suffix = ".tif";
print("\\Clear");
run("Close All");



////////////////////////// start template

input = getDirectory("Input directory");
// input = "C:\\Users\\sokip\\Dropbox (UMass Medical School)\\Microscopy\\Other\\test_Damon\\";
// output = getDirectory("Output directory");
output_path = input+File.separator+output;
File.makeDirectory(output_path);

/*
Dialog.create("File type");
Dialog.addString("File suffix: ", ".tif", 5);
Dialog.show();
suffix = Dialog.getString();
*/


print (macro_name);
print (output_path);
print ("Rearrange channels: " + rearrange);
print ("Minimum object size: " + min_size);
print ("C1 DOG (s1,s2), erode: " + C1_s1 + "," + C1_s2 + "," + C1_e);
print ("C2 DOG (s1,s2), erode: " + C2_s1 + "," + C2_s2 + "," + C2_e);
print ("Auto threshold: " + auto_thr);
print ("Min threshold: " + min_thr);
print ("Max threshold: " + max_thr);
print ("DAPI_gb: " + DAPI_gb);
print ("DAPI_dilate: " + DAPI_dilate);
print ("DAPI_thr: " + DAPI_thr);


print ("");
processFolder(input);
selectWindow("Log");
saveAs("text", output_path+File.separator+output);

function processFolder(input) {
	list = getFileList(input);
	count = 0;
	for (j=0; j<list.length; j++) {
		if(endsWith(list[j], suffix)) {
			print (list[j]);
			count++;
		}
	}
	print ("Found " + count + " files to process");
	print ("");

	if (dbug == true) {
		setBatchMode(false);
		/*
		if (count > 1) { setBatchMode(true); }
		else { setBatchMode(false); }
		*/
	} else { setBatchMode(true); }

	x = 1;
	for (i = 0; i < list.length; i++) {
		if ((i+1) % 5 == 0) { // run garbage collector every nth image
//			print ("GARBAGE COLLECTING");
			call("java.lang.System.gc");
		}
		// uncomment if you want to process folders recursively
		/*
		if(File.isDirectory(input + list[i]))
			processFolder("" + input + list[i]);
		*/
		if(endsWith(list[i], suffix)) {
			print ("Processing file " +x+ " out of " +count);
			processFile(input, output_path, list[i]);
			x++;
		}
	}
	
}

print ("DONE");
waitForUser("DONE!");


////////////////////////// end template
////////////////////////// start individual processing

function processFile(input, output_path, file) {
//	IJ.freeMemory(); // show available memory
	if (dbug == false) {
		run("Close All");
	}
	open(input+file);
	run("Select None");
/*
	if (rearrange != false) {
		run("Arrange Channels...", "new="+rearrange);
	}
*/	

	// Initiate variables
	base = getTitle();
	print (base);
	Stack.getDimensions(width, height, channels, slices, frames);

	// if image only has 2 channels, duplicate one of the channels
	if (channels == 2) {
		getVoxelSize(width, height, depth, unit);
		run("Make Substack...", "channels=2,2,1 slices=1-" + slices);
		setVoxelSize(width, height, depth, unit);
		channels++;
		close_title(base); // make substack makes a new image so we have to close one
		selectImage(1);
		rename(base);
	} else if (channels == 3) {
		if (rearrange != false) {
			run("Arrange Channels...", "new="+rearrange);
		}
	}
	
	run("Split Channels"); // expecting 3 channel image
	images = newArray(channels+1);
	proc = newArray(channels+1);
	titles = newArray(channels+1);
	nuc = newArray(channels);
	cyto = newArray(channels);
	nuc_obj = newArray(channels);
	cyto_obj = newArray(channels);
	
	
	// Parse separate channels into an array
	for (i=0; i<nImages(); i++) {
		selectImage(i+1);
		title = getTitle();
		print(title);
		if (startsWith(title,"C1")) { images[1] = getImageID(); } // skip images[0]
		else if (startsWith(title,"C2")) { images[2] = getImageID(); }
		else if (startsWith(title,"C3")) { images[3] = getImageID(); }
	}
	
	
	// Threshold RNAscope channels
	titles[1] = "" + DOG(images[1],C1_s1,C1_s2,C1_e);
	titles[2] = "" + DOG(images[2],C2_s1,C2_s2,C2_e);
	thresholds = newArray(3);
	thresholds[1] = calc_threshold(titles[1]);
	thresholds[2] = calc_threshold(titles[2]);


	// Threshold DAPI image
	selectImage(images[3]);
	run("Gaussian Blur...", "sigma="+DAPI_gb+" stack");
	run("8-bit");
	//setAutoThreshold("Otsu dark stack");
	//run("Convert to Mask", "method=Otsu background=Dark black"); // preserves 16-bit
	run("Auto Threshold", "method="+DAPI_thr+" ignore_black ignore_white white stack use_stack_histogram"); // ignore black and white pixels
	run("Fill Holes", "stack");

	// 2D dilate
	/*
	for (i=0; i<DAPI_dilate; i++) {
		run("Dilate", "stack");
	}
	*/

	// 3D dilate
	if (DAPI_dilate > 0) {
		selectImage(images[3]);
		temp_ID = getImageID();
		temp_title = getTitle();
		run("Morphological Filters (3D)", "operation=Dilation element=Ball x-radius=&DAPI_dilate y-radius=&DAPI_dilate z-radius=&DAPI_dilate");
		rename(temp_title);
		close_id(temp_ID);
		images[3] = getImageID();		
	}
	
	run("16-bit"); // Auto Threshold converts to 8-bit, so need to convert back to 16-bit
	run("Multiply...", "value=65535 stack");
	titles[3] = getTitle();
	
	
	// Close original images
	// close_id(images[1]);
	// close_id(images[2]);
	
	
	// Segment nuclear and cyotplasmic RNAscope signal
	nuc[1] = "" + nuc_cyto(titles[1],titles[3],"nuc");
	cyto[1] = "" + nuc_cyto(titles[1],titles[3],"cyto");

	nuc[2] = "" + nuc_cyto(titles[2],titles[3],"nuc");
	cyto[2] = "" + nuc_cyto(titles[2],titles[3],"cyto");

	
	// Count 3D objects
	print ("C1&C2 NUCLEAR AND CYOTPLASMIC OBJECTS:");
	summary = true;
	nuc_obj[1] = "" + obj_counter(nuc[1],thresholds[1],min_size,summary);
	cyto_obj[1] = "" + obj_counter(cyto[1],thresholds[1],min_size,summary);
	nuc_obj[2] = "" + obj_counter(nuc[2],thresholds[2],min_size,summary);
	cyto_obj[2] = "" + obj_counter(cyto[2],thresholds[2],min_size,summary);


/*
	//print ("C1&C2 NUCLEAR COLOC (OBJECT):");
	coloc_obj(nuc[1],nuc[2],min_size);
	
	//print ("C1&C2 CYTOPLASMIC COLOC (OBJECT):");
	coloc_obj(cyto[1],cyto[2],min_size);
*/


	if (save_files == true) {
//		save_dir = strip_ext(output+File.separator+base);
		save_dir = strip_ext(output_path+File.separator+base);
		File.makeDirectory(save_dir);
//		print ("save_dir: ", save_dir);
//		print ("base: ", base);
		n = nImages();
		for (i=1; i<=n; i++) {
			selectImage(i);
			title = getTitle();
			saveAs("tiff", save_dir+File.separator+title);
		}
	}

	if (save_stats == true) {
		for (i=1; i<=2; i++) {
			thr_here = thresholds[i];
			summary = false;
			// unsegmented image (nuc+cyto combined)
			selectWindow(titles[i]+".tif");
//			print ("SAVING STATS: " + titles[i]);
//			run("3D Objects Counter", "threshold=&thr_here slice=1 min.=&min_size max.=999999999 objects statistics");
			obj_counter(getTitle(), thr_here, min_size, summary);
			selectWindow("Statistics for " + titles[i] + ".tif");
			saveAs("Results", save_dir+File.separator+titles[i]+"_stats.tsv");

			// nuc image
			selectWindow(titles[i]+"_nuc.tif");
//			print ("SAVING STATS: " + titles[i]);
//			run("3D Objects Counter", "threshold=&thr_here slice=1 min.=&min_size max.=999999999 objects statistics");
			obj_counter(getTitle(), thr_here, min_size, summary);
			selectWindow("Statistics for " + titles[i] + "_nuc.tif");
			saveAs("Results", save_dir+File.separator+titles[i]+"_nuc_stats.tsv");

			// cyto image
			selectWindow(titles[i]+"_cyto.tif");
//			print ("SAVING STATS: " + titles[i]);
//			run("3D Objects Counter", "threshold=&thr_here slice=1 min.=&min_size max.=999999999 objects statistics");
			obj_counter(getTitle(), thr_here, min_size, summary);
			selectWindow("Statistics for " + titles[i] + "_cyto.tif");
			saveAs("Results", save_dir+File.separator+titles[i]+"_cyto_stats.tsv");
		}
	}
	print("");
	close_windows();
//	if (list.length > 1) { run("Close All"); }
}

////////////////// end individual processing






////////////////////////////////
//////////// FUNCTIONS
////////////////////////////////


// Difference of Gaussians function
function DOG (id, s1, s2, erode) {
	selectImage(id);
	title = getTitle();
	run("Duplicate...", "title=s"+s1+" duplicate");
	s1_title = getTitle();
	run("Duplicate...", "title=s"+s2+" duplicate");
	s2_title = getTitle();
	
	selectWindow(s1_title);
	run("Gaussian Blur...", "sigma="+s1+" stack");
	selectWindow(s2_title);
	run("Gaussian Blur...", "sigma="+s2+" stack");
	
	imageCalculator("Subtract create stack", s1_title, s2_title);
/*
//	print(getTitle());
//	run("8-bit");
//	run("Watershed", "stack");
//	run("16-bit");
	run("Multiply...", "value=65535 stack");	
//	temp = getImageID();
*/

//	selectImage(temp);
	
	for (i=0; i<erode; i++) {
		run("8-bit");
		run("Erode", "stack");
//		print("ERODE");
	}
	if (endsWith(title,".tif")) {
		title = substring(title,0,lastIndexOf(title,".tif"))+"_"+s1_title+"-"+s2_title;
	} else {
		title += "_"+s1_title+"-"+s2_title;
	}
	title += "_e"+erode;
	rename(title);
	close_title(s1_title);
	close_title(s2_title);
	return title;
}


// Nucleus or cytoplasm segmenter
function nuc_cyto(image, reference, area) {
	selectWindow(image);
	title = getTitle() + "_" + area;
	if (area == "nuc") { operation = "AND"; }
	if (area == "cyto") { operation = "Subtract"; }

	imageCalculator(operation + " create stack", image, reference);
	rename(title);
	return title;
}


// Measure colocalization on PIXEL basis
function coloc_pixels (image1, image2) {
	imageCalculator("AND create stack", image1, image2);
	return (measure_stack_intensity());
}

/*
// Measure colocalization on OBJECT basis
function coloc_obj (image1, image2, min_size) {
	imageCalculator("AND create stack", image1, image2);
	title = getTitle();
	obj_counter(title, min_size);
	return title;
}
*/

// Measure stack intensity of current image
function measure_stack_intensity () {
	intensity = 0;
	setOption("Stack position", true);
	for (n=1; n<=nSlices; n++) {
		setSlice(n);
		run("Measure");
		intensity += getResult("RawIntDen");
	}
	selectWindow("Results");
	run("Close");
	return intensity;
}

/*
// 3D objects counter macro
function obj_counter (title,thr,min_size,summary) {
	selectWindow(title);
	Stack.getStatistics(voxelCount, mean, min, max, stdDev);
	if (max < thr) { // don't run 3DOC if the threshold is above the maximum pixel value in the image
		print (getTitle() + ": 0 objects detected (Size filter set to " + min_size + " voxels, threshold set to: " + thr + ")");
	} else if (summary == true) {
		run("3D Objects Counter", "threshold=&thr slice=1 min.=&min_size max.=999999999 objects statistics summary");	
	} else if (summary == false) {
		run("3D Objects Counter", "threshold=&thr slice=1 min.=&min_size max.=999999999 objects statistics");
	}
//	setThreshold(thr, 65535);
//	run("Subtract...", "value=&thr stack");
//	run("Multiply...", "value=65535 stack");
	return title;
}
*/



// 3D objects counter macro
function obj_counter (title,thr,min_size,summary) {
	selectWindow(title);
	Stack.getStatistics(voxelCount, mean, min, max, stdDev);
	if (max < thr) { thr = max; }
	if (summary == true) {
		run("3D Objects Counter", "threshold=&thr slice=1 min.=&min_size max.=999999999 objects statistics summary");
	} else if (summary == false) {
		run("3D Objects Counter", "threshold=&thr slice=1 min.=&min_size max.=999999999 objects statistics");
	}
//	setThreshold(thr, 65535);
//	run("Subtract...", "value=&thr stack");
//	run("Multiply...", "value=65535 stack");
	return title;
}


function calc_threshold (title) {
	selectWindow(title);
	setAutoThreshold(auto_thr + " dark no-reset stack");
	getThreshold(thr, upper);
	if (thr < min_thr) {
		thr = min_thr;
	} else if (thr > max_thr) {
		thr = max_thr;
	}
	return thr;
}


// Close image by id
function close_id(id) {
	selectImage(id);
	close();
}

// Close image by title
function close_title(title) {
	selectWindow(title);
	close();
}

// Close all non-image windows EXCEPT "Log"
function close_windows() {
	list = getList("window.titles");
	for (i=0; i<list.length; i++) {
		winame = list[i];
		selectWindow(winame);
		if (winame != "Log") { run("Close"); }
	}
}


// Strip extension from a string by looking for last "." character
function strip_ext (s) {
	if (indexOf(s,".") > 0) {
		return substring(s,0,lastIndexOf(s,"."));
	}
	return s;
}