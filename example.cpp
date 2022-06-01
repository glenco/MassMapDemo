#include <slsimlib.h>
#include <sstream>
#include <iomanip>
#include <omp.h>
#include <thread>
#include <mutex>
#include "grid_maintenance.h"
#include "gridmap.h"


int main(int arg,char **argv){

  /********************
   set parameter file name:
   A default name is used or a name is taken as a command line argument
   *********************/

//  std::string paramfile;
//  std::cout << "initializing model" << std::endl;
//  //string paramfile;
//  if(arg > 1) paramfile.assign(argv[1],strlen(argv[1]));
//  else paramfile = "sample_paramfile";
//  std::cout << "using parameter file: " << paramfile << std::endl;
//
  // read parameter file
  //InputParams params(paramfile);

  long seed = -1827674;
  COSMOLOGY cosmo(CosmoParamSet::Planck18);
  std::cout << "Lens constructed" << std::endl;
  
  //double center[] = {0.3*pi/180,-0.25*pi/180};
  double center[] = {0,0};  // center of grids
  double range = 0.3*PI/180; // range of grids in radians
  size_t Ninit = 1024; // the initial number of pixels to a side in the grid
  double zsource = 2,zlens = 0.5;
  
  /********************
   construct lens:
   *******************/

  // make an empty Lens
  Lens lens(&seed,zsource,cosmo);
  {
    // make a mass map LensHalo
    double mass_units = 14204545454.5455/1.13572841655696E-05;
    //LensHaloMassMap map("map.fits",,zlens,mass_units,1,false,cosmo);
    LensHaloMassMap map("map.fits",PixelMapType::pix_map, 1, true, cosmo);
    
    double Sigma_crit = cosmo.SigmaCrit(zlens,zsource);
    double Dl = cosmo.angDist(zlens);
    double resolution_in_mpc = map.getRangeMpc()/map.getNx();
    
    // Print an image of the surface density of the LensHalo
    PixelMap image = map.map_variables(LensingVariable::KAPPA,1000,1000,resolution_in_mpc);
    image *= 1/Sigma_crit;  // rescale to convergence
    image.printFITS("!kappa_image.fits");
    
    // Print an image of the first component of the LensHalo
    image = map.map_variables(LensingVariable::ALPHA1,1000,1000,resolution_in_mpc);
    image *= 1/Sigma_crit/Dl; // rescale to radians
    image.printFITS("!alpha1_image.fits");
    
    // Print an image the photenital of the LensHalo
    image = map.map_variables(LensingVariable::PHI,1000,1000,resolution_in_mpc);
    image *= 1/Sigma_crit/Dl/Dl; // rescale to radians^2
    image.printFITS("!phi_image.fits");

    // These output images are fine for a single plane lens
    // but will not be valid for a multiplane lens.
    
    lens.moveinMainHalo(map, true);
  }

  /** 
   Here a uniform grid is constructed that is not capable 
   of refinement for images, caustics, etc.  This is the 
   simplest way to make a map.
   In the construction the rays are shot through the lens 
   in parallel if the N_THREAD is set to multiple threads 
   in the GLAMER build.
  **/
  {
    std::cout << "Constructing initial GridMap ..." << std::endl;
    GridMap gridmap(&lens,Ninit,center,range);
    std::cout << "constructed" << std::endl;

    /***
     make an images of kappa
     The center,range and pixel numbers do not have to match the grid, but 
     in this case they are set to match.
     ****/
    gridmap.writeFitsUniform(center,gridmap.getInitNgrid()
                           ,gridmap.getInitNgrid(),LensingVariable::KAPPA,"!gridmap");

  }
  /****************************
   Here a Grid is constucted which can be refined for finding images, etc.
   ****************************/
  
  std::cout << "Constructing initial Grid ..." << std::endl;
  Grid grid(&lens,Ninit,center,range);
  std::cout << "constructed" << std::endl;
  
  /*
   There we make same maps of the lensing quantities.
   This could be done with GridMap since no refinement has been done yet.
   The "!" infront of the file name causes it to overwrite a file with that name.  Suffixes are added (eg .kappa.fits).
   */
  grid.writeFits(center,grid.getInitNgrid(), grid.getInitRange()/grid.getInitNgrid()
                 ,LensingVariable::KAPPA,"!initgrid");
  grid.writeFits(center,grid.getInitNgrid(), grid.getInitRange()/grid.getInitNgrid()
                 ,LensingVariable::ALPHA,"!initgrid");
  grid.writeFits(center,grid.getInitNgrid(), grid.getInitRange()/grid.getInitNgrid()
                 ,LensingVariable::ALPHA1,"!initgrid");
  grid.writeFits(center,grid.getInitNgrid(), grid.getInitRange()/grid.getInitNgrid()
                 ,LensingVariable::ALPHA2,"!initgrid");
  grid.writeFits(center,grid.getInitNgrid(), grid.getInitRange()/grid.getInitNgrid()
                 ,LensingVariable::GAMMA,"!initgrid");
  
  /******************************************
   Now we are going to look for same caustics
   ******************************************/

  // The CriticalCurve class contains information about a caustic
  std::vector<ImageFinding::CriticalCurve> critcurves(100);
  int Ncrit;  // number of caustics
  // resolution to which the critical curves should be refined (radians)
  PosType resolution = 0.05*PI/180/60/60;
  
  std::cout << "Looking for critical curves ..." << std::endl;
  // Find all the critical curves
  ImageFinding::find_crit(&lens,&grid,critcurves, &Ncrit,resolution);
  std::cout << Ncrit << " critical curves found." << std::endl;

  PosType Xrange[2]={0,0},Yrange[2]={0,0};
  
  if(Ncrit > 0){
    Point_2d p1,p2;
    
    // find a box on the image plane that cantains all of the critical curves
    for(int ii = 0;ii< Ncrit;++ii){
      critcurves[ii].CritRange(p1,p2);
 
      Xrange[0] = MIN(Xrange[0],p1[0]);
      Xrange[1] = MAX(Xrange[1],p2[0]);
      
      Yrange[0] = MIN(Yrange[0],p1[1]);
      Yrange[1] = MAX(Yrange[1],p2[1]);
    }
    
    // make a PixelMap which is used for IO of images
    PixelMap map(center
                 ,(size_t)(MAX(Xrange[1]-Xrange[0],Yrange[1]-Yrange[0])/resolution/3)
                 ,resolution*3);
    
    // this draws the critical curves on the image
    for(int ii = 0;ii< Ncrit;++ii){
      map.AddCurve(critcurves[ii].critcurve,ii+1);
    }
    
    /**************************
     Now we are going to pick a caustic, put a source 
     in it and find its images
     **************************/
    int nc=0;  // pick a caustic
    // draw it on the image
    map.AddCurve(critcurves[nc].caustic_curve_outline,2);

    std::vector<Point_2d> ys;
    // this is a random number generator
    Utilities::RandomNumbers_NR rng(seed);
    
    // This finds random points within this caustic.
    // In this case it is just one point.
    critcurves[nc].RandomSourceWithinCaustic(1,ys,rng);

    // ImageInfo is a class that contains information about images
    std::vector<ImageInfo> imageinfo(100);
    int Nimages;  // number of images that will be found
    size_t Nimagepoints;  // total number of points within the images
    
    // this finds the images made by a source at ys[0] and stores
    // information about them in imageinfo.
    // These are just circular sources.
    //  More realistic sources can be mapped by construction a source
    //  and using ImageFinding::map_imagesISOP().
    ImageFinding::find_images_kist(&lens,ys[0].x,1.1*PI/180/60/60
                    ,&grid,&Nimages,imageinfo,&Nimagepoints,
                            0, true, 2);
    
    std::cout << "Number of images: " << Nimages << std::endl;
    // add images to the PixelMap
    map.AddImages(imageinfo.data(),Nimages,-1);
    
    // ooutput the PixelMap as a fits file
    map.printFITS("!test.fits");
  }
  
  // If the grid is now output at twice the original resolution
  // some additional structure might be see because of the refinement
  grid.writeFits(center,2*grid.getInitNgrid(),grid.getInitRange()/grid.getInitNgrid()/2,
                 LensingVariable::INVMAG,"!initgrid_refined");
  grid.writeFits(center,2*grid.getInitNgrid(),grid.getInitRange()/grid.getInitNgrid()/2, LensingVariable::KAPPA,"!initgrid_refined");
  
  
  /*************************************************************
   Many other things are possible and easily done with GLAMER.  
   Read the documantation for a more complete sicription of functionality.
   *************************************************************/
 }
