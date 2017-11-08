path= newArray(16);
//path[0]="/Volumes/Birdbox III/Aptamer_studies/13_02 Cells/20170213_AptPaint_Cell_GS/561_1/"
path[0]="D:/Dropbox (Cambridge University)/Group/Livesey_LysateIII_2017-08-21_16-17-26/H9-T100/561/";
path[1]="D:/Dropbox (Cambridge University)/Group/Livesey_LysateIII_2017-08-21_16-17-26/H9-T300/561/";
path[2]="D:/Dropbox (Cambridge University)/Group/Livesey_LysateIII_2017-08-21_16-17-26/H9-V100/561/";
path[3]="D:/Dropbox (Cambridge University)/Group/Livesey_LysateIII_2017-08-21_16-17-26/H9-V300/561/";
path[4]="D:/Dropbox (Cambridge University)/Group/Livesey_LysateIII_2017-08-21_16-17-26/I4-T100/561/";
path[5]="D:/Dropbox (Cambridge University)/Group/Livesey_LysateIII_2017-08-21_16-17-26/I4-T300/561/";
path[6]="D:/Dropbox (Cambridge University)/Group/Livesey_LysateIII_2017-08-21_16-17-26/I4-V100/561/";
path[7]="D:/Dropbox (Cambridge University)/Group/Livesey_LysateIII_2017-08-21_16-17-26/I4-V300/561/";
path[8]="D:/Dropbox (Cambridge University)/Group/Livesey_LysateIII_2017-08-21_16-17-26/M146I-T100/561/";
path[9]="D:/Dropbox (Cambridge University)/Group/Livesey_LysateIII_2017-08-21_16-17-26/M146I-T300/561/";
path[10]="D:/Dropbox (Cambridge University)/Group/Livesey_LysateIII_2017-08-21_16-17-26/M146I-V100/561/";
path[11]="D:/Dropbox (Cambridge University)/Group/Livesey_LysateIII_2017-08-21_16-17-26/M146I-V300/561/";
path[12]="D:/Dropbox (Cambridge University)/Group/Livesey_LysateIII_2017-08-21_16-17-26/ND-T100/561/";
path[13]="D:/Dropbox (Cambridge University)/Group/Livesey_LysateIII_2017-08-21_16-17-26/ND-T300/561/";
path[14]="D:/Dropbox (Cambridge University)/Group/Livesey_LysateIII_2017-08-21_16-17-26/ND-V100/561/";
path[15]="D:/Dropbox (Cambridge University)/Group/Livesey_LysateIII_2017-08-21_16-17-26/ND-V300/561/";

for (i=0; i<path.length; i++){
k=0;
count = 0;	
n=0;
dir=path[i];
processFiles(dir);

}


/////////// Function to process the files ///////////////

function processFiles(dir) 

	{

	list = getFileList(dir);

		for (i=0; i<list.length; i++)

		 {

		if (!startsWith(list[i],"Log"))

			{

			if (endsWith(list[i], "/"))

		              processFiles(""+dir+list[i]);

         			 else 

			{

		             showProgress(n++, count);

            			path = dir+list[i];

	            		processFile(path);

			}

			}

		}

	}

function processFile(path) 

	{
		       	if (endsWith(path, ".tif") || endsWith(path, ".tiff") && !endsWith(path, "SuperRes.tif")) 
 
		{
			//open(path);

			run("TIFF Virtual Stack...", "open=["+path+"]");

			
			name = getTitle(); 
            dotIndex = indexOf(name, "."); 
            file= substring(name, 0, dotIndex); 
			
			rename("Image");
			myDir = dir+k+"/";
			print(myDir);

  			File.makeDirectory(myDir);

			run("Peak Fit", "config_file=[C:\\Users\\DRW\\gdsc.smlm.settings.xml] calibration=131.5 gain=55.50 exposure_time=25 initial_stddev0=2.000 initial_stddev1=2.000 initial_angle=0.000 smoothing=0.50 smoothing2=3 search_width=3 fit_solver=[Least Squares Estimator (LSE)] fit_function=Circular fit_criteria=[Least-squared error] significant_digits=5 coord_delta=0.0001 lambda=10.0000 max_iterations=20 fail_limit=10 include_neighbours neighbour_height=0.30 residuals_threshold=1 duplicate_distance=0.50 shift_factor=2 signal_strength=30 width_factor=2 precision=20 results_table=Uncalibrated image=[Localisations (width=precision)] weighted equalised image_precision=5 image_scale=8 results_dir=["+myDir+"] stack");
			
			selectWindow("Image (LSE) SuperRes");
	//		saveAs("Tiff", ""+myDir+"SuperRes.tif");
			close();
			close();
			selectWindow("Log");
 			selectWindow("Fit Results");
			saveAs("text",myDir+"FitResults.txt");
	k++;
 			
run("Close"); 

		}
	}