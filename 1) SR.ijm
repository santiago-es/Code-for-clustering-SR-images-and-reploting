path= newArray(13);
//path[0]="/Volumes/Birdbox III/Aptamer_studies/13_02 Cells/20170213_AptPaint_Cell_GS/561_1/"

path[0]="/Volumes/BirdBox2/20170414_Cell_exps/AB_axons/4/561/"
path[1]="/Volumes/BirdBox2/20170414_Cell_exps/Media_2017-04-13_18-15-33/AB_axoncontrol/561/"
path[2]="/Volumes/BirdBox2/20170414_Cell_exps/Media_2017-04-13_18-15-33/AB_body/561"
path[3]="/Volumes/BirdBox2/20170414_Cell_exps/Media_2017-04-13_18-15-33/AB_bodycontrol/"
path[4]="/Volumes/BirdBox2/20170414_Cell_exps/AB_bodies/1/561/"
path[5]="/Volumes/BirdBox2/20170414_Cell_exps/AB_bodies/2/561/"
path[6]="/Volumes/BirdBox2/20170414_Cell_exps/AB_bodies/3/561/"
path[7]="/Volumes/BirdBox2/20170414_Cell_exps/AB_bodies/4/561/"
path[8]="/Volumes/BirdBox2/20170414_Cell_exps/AB_bodies/5/561/"
path[9]="/Volumes/BirdBox2/20170414_Cell_exps/AB_axons/1/561/"
path[10]="/Volumes/BirdBox2/20170414_Cell_exps/AB_axons/2/561/"
path[11]="/Volumes/BirdBox2/20170414_Cell_exps/AB_axons/3/561/"
path[12]="/Volumes/BirdBox2/20170414_Cell_exps/Media_2017-04-13_18-15-33/AB_axon/561/"

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

			run("Peak Fit", "config_file=[/Users/Mathew/gdsc.smlm.settings.xml] calibration=131.5 gain=55.50 exposure_time=25 initial_stddev0=2.000 initial_stddev1=2.000 initial_angle=0.000 smoothing=0.50 smoothing2=3 search_width=3 fit_solver=[Least Squares Estimator (LSE)] fit_function=Circular fit_criteria=[Least-squared error] significant_digits=5 coord_delta=0.0001 lambda=10.0000 max_iterations=20 fail_limit=10 include_neighbours neighbour_height=0.30 residuals_threshold=1 duplicate_distance=0.50 shift_factor=2 signal_strength=30 width_factor=2 precision=20 results_table=Uncalibrated image=[Localisations (width=precision)] weighted equalised image_precision=5 image_scale=8 results_dir=["+myDir+"] stack");
			
			selectWindow("Image (LSE) SuperRes");
			saveAs("Tiff", ""+myDir+"SuperRes.tif");
			close();
			close();
			selectWindow("Log");
 			selectWindow("Fit Results");
			saveAs("text",myDir+"FitResults.txt");
	k++;
 			
run("Close"); 

		}
	}