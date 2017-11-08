filelist= newArray(9);
filelist[0]="/Volumes/Birdbox III/Aptamer_studies/Timecourse_PAINT/0h/561/0/"
filelist[1]="/Volumes/Birdbox III/Aptamer_studies/Timecourse_PAINT/0h/561/1/"
filelist[2]="/Volumes/Birdbox III/Aptamer_studies/Timecourse_PAINT/2h/apt/0/"
filelist[3]="/Volumes/Birdbox III/Aptamer_studies/Timecourse_PAINT/2h/apt/1/"
filelist[4]="/Volumes/Birdbox III/Aptamer_studies/Timecourse_PAINT/4h/apt/0/"
filelist[5]="/Volumes/Birdbox III/Aptamer_studies/Timecourse_PAINT/4h/apt/1/"
filelist[6]="/Volumes/Birdbox III/Aptamer_studies/Timecourse_PAINT/6h/apt/0/"
filelist[7]="/Volumes/Birdbox III/Aptamer_studies/Timecourse_PAINT/6h/apt/1/"
filelist[8]="/Volumes/Birdbox III/Aptamer_studies/Timecourse_PAINT/8h/apt/0/"
filelist[9]="/Volumes/Birdbox III/Aptamer_studies/Timecourse_PAINT/8h/apt/1/"
filelist[10]="/Volumes/Birdbox III/Aptamer_studies/Timecourse_PAINT/22h/apt/0/"
filelist[11]="/Volumes/Birdbox III/Aptamer_studies/Timecourse_PAINT/22h/apt/1/"

for (i=0; i<filelist.length; i++){


run("Text Image... ", "open=["+filelist[i]+"imagegauss2.txt]");
run("Rotate 90 Degrees Right");
run("Flip Horizontally");

run("Red Hot");
setMinAndMax(0, 30);
saveAs("Tiff", ""+filelist[i]+"SR_clustersonly.tif");


close();

}



