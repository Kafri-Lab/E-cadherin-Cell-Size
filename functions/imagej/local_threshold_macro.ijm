
//input = getDirectory("Input directory");
//output = getDirectory("Output directory");
input = "Z:/DanielS/Images/zoo_animal/hepatocyte_images/old/ALL_backup/"
output = "Z:/DanielS/Images/zoo_animal/hepatocyte_images/old/Phansalkar_threshold_nuc/"

processFolder(input);
 
function processFolder(input) {
    list = getFileList(input);
    for (i = 0; i < list.length; i++) {
        processFile(input, output, list[i], "Phansalkar");
    }
    print("Done");
}
 
function processFile(input, output, file, thresh) {
    // do the processing here by replacing
    // the following two lines by your own code
    print("Processing: " + input + file);
    open(input + file);
	run("Split Channels");
	selectWindow(file + " (blue)");
	run("Median...", "radius=9");
	run("Options...", "iterations=1 count=1 black do=Nothing");
    run("Auto Local Threshold", "method=Phansalkar radius=30 parameter_1=0 parameter_2=0 white");
    run("Fill Holes");
    saveAs("TIFF", output+file+"_thresh");
    close();
}
