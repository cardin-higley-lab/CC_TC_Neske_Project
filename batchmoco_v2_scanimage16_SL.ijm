/////batch motion corrects series of .tif images in input directory with moco, and outputs results to output directory/////

//user-defined parameters for input/output directories and for moco and reference
input=getDirectory("Choose a directory with images to be motion corrected");   //directory with images to be motion corrected
output=getDirectory("Choose an output directory");   //directory for output motion corrected images

output_details=output+"moco_info/";   //directory for output .xls translation coord and ref frame
if (!File.exists(output_details)) File.makeDirectory(output_details); 

//initial guesses for reference frame details and moco parameter
first_ref_frame_init=1; 
last_ref_frame_init=200;
moco_max_dist_init=20;

//check if have reference frame to use. If not, get parameters to create one.  Get moco param
hasRef=getBoolean("Do you have a saved reference image that you would like to use?");
if (hasRef==0){
   Dialog.create("Choose parameters");
   Dialog.addString("Index of first reference frame:",first_ref_frame_init); 
   Dialog.addString("Index of last reference frame:",last_ref_frame_init); 
   Dialog.addString("The maximum distance (in pixels) to be translated in the x and y directions:",moco_max_dist_init);
   Dialog.show();
   first_ref_frame=Dialog.getString();
   last_ref_frame=Dialog.getString();
   moco_max_dist=Dialog.getString();
} else {
   waitForUser( "Pause","Open your reference frame and press OK when ready"); 
   ref_frame=getTitle;
   saveAs("Tiff", output_details+ref_frame);
   selectWindow(ref_frame); close();
   moco_max_dist=getNumber("The maximimum distance (in pixels) to be translated in the x and y directions",moco_max_dist_init);    // max distance to use for moco
}

/////////////////////////////////////////////

//if not given a reference frame, make ref frame from first video.  If does, save. Moco all files in input dir

list=getFileList(input);
tiffiles=newArray(0);
for (listind=0; listind<list.length; listind++){
     if(endsWith(list[listind],".tif")) tiffiles=Array.concat(tiffiles,list[listind]);
}
for (fileind=0; fileind<tiffiles.length; fileind++){
    if (fileind==0 && hasRef==0) {
         ref_frame=make_ref_frame(tiffiles[fileind]);   
    }
    run_batch_moco(tiffiles[fileind],ref_frame);
}

////////////////////////////////// 

//makes reference frame using average intensity over a defined number of ref frames.  Names ref image, saves and closes
function make_ref_frame(filename){
    run("Bio-Formats", "open=[" + input + filename + "] autoscale color_mode=Default view=Hyperstack stack_order=XYCZT");
    run("Z Project...", "start=first_ref_frame stop=last_ref_frame projection=[Average Intensity]");
    filename_noext=dropext(filename);
    ref_frame= filename_noext + "_REF_avg" +first_ref_frame+ "-" + last_ref_frame+ ".tif";     //name ref frame
    saveAs("Tiff", output_details + ref_frame);
    selectWindow(ref_frame); close();
    selectWindow(filename); close();
    return ref_frame;
}

//runs moco using reference frame for first video
function run_batch_moco(filename,ref_frame){
    run("Bio-Formats", "open=[" + input + filename + "] autoscale color_mode=Default view=Hyperstack stack_order=XYCZT");
    run("Bio-Formats", "open=[" + output_details + ref_frame + "] autoscale color_mode=Default view=Hyperstack stack_order=XYCZT");
    ref_frame1= "moco_info/"+ref_frame;     //name ref frame (SL Added)
    run("moco ", "value=moco_max_dist downsample_value=0 template="+ref_frame1+" stack="+filename+" log=[Generate log file] plot=[No plot]"); // SL modified 
    selectWindow(filename); close();  
    run("Bio-Formats", "open=[" + input + filename + "] autoscale color_mode=Default view=Hyperstack stack_order=XYCZT");
   translate_xy();
   save_and_close(filename,ref_frame);
}

/////////////////////////////////////////////////////////////
// some helper functions

function translate_xy(){
    for (i=1;i<=nSlices;i++){
        xtrans=getResult("x",i-1);
        ytrans=getResult("y",i-1);
        setSlice(i);
        run("Translate...", "x="+xtrans+" y="+ytrans+" interpolation=None slice"); 
    }
}

//defines names for motion corrected images and translation coordinates (in .xls file), saves and closes 
function save_and_close(filename,ref_frame){
    filename_noext=dropext(filename);
    if (hasRef==0){
        moco_name=filename_noext+"_moco"+moco_max_dist+"_ref"+first_ref_frame+ "-" +last_ref_frame;
        xytrans_name=filename_noext+"_xytrans_moco"+moco_max_dist+"_ref"+first_ref_frame+ "-" +last_ref_frame+ ".xls";
    } else {
        moco_name=filename_noext+"_moco"+moco_max_dist;
        xytrans_name=filename_noext+"_xytrans_moco"+moco_max_dist+".xls";
    }
    selectWindow("New Stack"); close();
    selectWindow(filename); saveAs("Tiff", output + moco_name); close();
    selectWindow(ref_frame1); close();
    selectWindow("Results"); saveAs("Results", output_details+xytrans_name); run("Close");
}


function dropext(filename){
    dotIndex = indexOf(filename, ".");
    filename_noext = substring(filename, 0, dotIndex); 
    return filename_noext;
} 
