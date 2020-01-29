// Function that computes the 1D cubic interpolation. g0,g1,g2,g3 contain the values of the voxels that we want to interpolate
// factor contains the distance from the value g1 to the position we want to interpolate to.
// We assume the out of bounce checks have been done earlier
public float cubicinterpolate(float g0,float g1,float g2,float g3,float factor){
        float result=0.0f;

        // computing according to the formula in the slide
        // define it with parameters: factor and a
        // a was chosen to be -0.75
        float[]coeffs={0.0f,0.0f,0.0f,0.0f};
        coeffs[0]=((a*(factor+1)-5*a)*(factor+1)+8*a)*(factor+1)-4*a;
        coeffs[1]=((a+2)*factor-(a+3))*factor*factor+1;
        coeffs[2]=((a+2)*(1-factor)-(a+3))*(1-factor)*(1-factor)+1;
        coeffs[3]=1.f-coeffs[0]-coeffs[1]-coeffs[2];

        return result=(float)(coeffs[0]*g0+coeffs[1]*g1+coeffs[2]*g2+coeffs[3]*g3);

}


// 2D cubic interpolation implemented here. We do it for plane XY. Coord contains the position.
// We assume the out of bounce checks have been done earlier
public float bicubicinterpolateXY(double[]coord,int z){

        // get the coordinations
        int x=(int)Math.floor(coord[0]);
        int y=(int)Math.floor(coord[1]);

        // compute the factor to be computed in the cubicinterpolate function
        float fac_x=(float)coord[0]-x;
        float fac_y=(float)coord[1]-y;

        // keep z as the same
        // operate the cubicinterpolation for plane XY
        float t0=cubicinterpolate(getVoxel(x-1,y-1,z),getVoxel(x-1,y,z),getVoxel(x-1,y+1,z),getVoxel(x-1,y+2,z),fac_y);
        float t1=cubicinterpolate(getVoxel(x,y-1,z),getVoxel(x,y,z),getVoxel(x,y+1,z),getVoxel(x,y+2,z),fac_y);
        float t2=cubicinterpolate(getVoxel(x+1,y-1,z),getVoxel(x+1,y,z),getVoxel(x+1,y+1,z),getVoxel(x+1,y+2,z),fac_y);
        float t3=cubicinterpolate(getVoxel(x+2,y-1,z),getVoxel(x+2,y,z),getVoxel(x+2,y+1,z),getVoxel(x+2,y+2,z),fac_y);

        float result=cubicinterpolate(t0,t1,t2,t3,fac_x);

        return result;

}

// 3D cubic interpolation implemented here given a position in the volume given by coord.
public float getVoxelTriCubicInterpolate(double[]coord){
        if(coord[0]< 1||coord[0]>(dimX-3)||coord[1]< 1||coord[1]>(dimY-3)
            ||coord[2]< 1||coord[2]>(dimZ-3)){
            return 0;
        }

        // define the factor to be used in the final dimension
        int z=(int)Math.floor(coord[2]);
        float fac_z=(float)coord[2]-z;

        // add z axis and make it to be 3D cubic interpolation
        // the result is based on 2D cubic interpolation
        float t0=bicubicinterpolateXY(coord,z-1);
        float t1=bicubicinterpolateXY(coord,z);
        float t2=bicubicinterpolateXY(coord,z+1);
        float t3=bicubicinterpolateXY(coord,z+2);

        float result=cubicinterpolate(t0,t1,t2,t3,fac_z);
        if(result<0)
            result=0;

        return result;

}


//Function that updates the "image" attribute (result of renderings)
// using the slicing technique.
public void slicer(double[]viewMatrix){

        // we start by clearing the image
        resetImage();

        // vector uVec and vVec define the view plane,
        // perpendicular to the view vector viewVec which is going from the view point towards the object
        // uVec contains the up vector of the camera in world coordinates (image vertical)
        // vVec contains the horizontal vector in world coordinates (image horizontal)
        double[]viewVec=new double[3];
        double[]uVec=new double[3];
        double[]vVec=new double[3];
        getViewPlaneVectors(viewMatrix,viewVec,uVec,vVec);

        // The result of the visualization is saved in an image(texture)
        // we update the vector according to the resolution factor
        // If the resolution is 0.25 we will sample 4 times more points.
        for(int k=0;k< 3;k++){
            uVec[k]=res_factor*uVec[k];
            vVec[k]=res_factor*vVec[k];
        }

        // compute the volume center
        double[]volumeCenter=new double[3];
        computeVolumeCenter(volumeCenter);

        // Here will be stored the 3D coordinates of every pixel in the plane
        double[]pixelCoord=new double[3];

        // We get the size of the image/texture we will be puting the result of the
        // volume rendering operation.
        int imageW=image.getWidth();
        int imageH=image.getHeight();

        int[]imageCenter=new int[2];
        // Center of the image/texture
        imageCenter[0]=imageW/2;
        imageCenter[1]=imageH/2;

        // imageW/ image H contains the real width of the image we will use given the resolution.
        //The resolution is generated once based on the maximum resolution.
        imageW=(int)(imageW*((max_res_factor/res_factor)));
        imageH=(int)(imageH*((max_res_factor/res_factor)));

        // sample on a plane through the origin of the volume data
        double max=volume.getMaximum();

        // Color that will be used as a result
        TFColor pixelColor=new TFColor();
        // Auxiliar color
        TFColor colorAux;

        // Contains the voxel value of interest
        int val;

        //Iterate on every pixel
        for(int j=imageCenter[1]-imageH/2;j<imageCenter[1]+imageH/2;j++){
            for(int i=imageCenter[0]-imageW/2;i<imageCenter[0]+imageW/2;i++){

                // computes the pixelCoord which contains the 3D coordinates of the pixels (i,j)
                computePixelCoordinatesFloat(pixelCoord,volumeCenter,uVec,vVec,i,j);

                //we now have to get the value for the in the 3D volume for the pixel
                //we can use a nearest neighbor implementation like this:
                //val = volume.getVoxelNN(pixelCoord);
                //you have also the function getVoxelLinearInterpolated in Volume.java
                //val = (int) volume.getVoxelLinearInterpolate(pixelCoord);
                //you have to implement this function below to get the cubic interpolation
                val=(int)volume.getVoxelTriCubicInterpolate(pixelCoord);

                // Map the intensity to a grey value by linear scaling
                pixelColor.r=(val/max);
                pixelColor.g=pixelColor.r;
                pixelColor.b=pixelColor.r;

                // the following instruction makes intensity 0 completely transparent and the rest opaque
                // pixelColor.a = val > 0 ? 1.0 : 0.0;   compo
                // Alternatively, apply the transfer function to obtain a color using the tFunc attribute
                // colorAux= tFunc.getColor(val);
                // pixelColor.r=colorAux.r;pixelColor.g=colorAux.g;pixelColor.b=colorAux.b;pixelColor.a=colorAux.a;
                // IMPORTANT: You can also simply use pixelColor = tFunc.getColor(val); However then you copy by reference and this means that if you change
                // pixelColor you will be actually changing the transfer function So BE CAREFUL when you do this kind of assignments
                //BufferedImage/image/texture expects a pixel color packed as ARGB in an int
                //use the function computeImageColor to convert your double color in the range 0-1 to the format need by the image
                int pixelColor_i=computeImageColor(pixelColor.r,pixelColor.g,pixelColor.b,pixelColor.a);
                image.setRGB(i,j,pixelColor_i);
                }
            }
}


// Compute the opacity based on the value of the pixel and the values of the
// triangle widget tFunc2D contains the values of the baseintensity and radius
// tFunc2D.baseIntensity, tFunc2D.radius they are in image intensity units
public double computeOpacity2DTF(double material_value,double material_r,double voxelValue,double gradMagnitude){

        double opacity=0.0;

        // get basic opacity from tFunc2D
        TFColor baseColor=tFunc2D.color;

        // define components to be used
        // compute the material angle and voxel angle
        double maxGrad=gradients.getMaxGradientMagnitude();
        double angle=Math.atan(material_r/maxGrad);
        double opposite=Math.abs(voxelValue-material_value);
        double adjacent=gradMagnitude;
        double voxelAngle=Math.atan(opposite/adjacent);

        // if the material equals the voxel and the gradient is zero, use basic opacity
        // else if the voxel angle is smaller, use another strategy
        if(material_value==voxelValue&&gradMagnitude==0){
            opacity=baseColor.a;
        }else if(voxelAngle<angle){
            opacity=baseColor.a*((voxelAngle/angle));
        }else{
            opacity=0.0;
        }

        return opacity;
}


//Function that updates the "image" attribute using the compositing// accumulatted raycasting
//It returns the color assigned to a ray/pixel given it's starting point (entryPoint) and the direction of the ray(rayVector).
// exitPoint is the last point.
//ray must be sampled with a distance defined by the sampleStep
public int traceRayComposite(double[]entryPoint,double[]exitPoint,double[]rayVector,double sampleStep){
        double[]lightVector=new double[3];

        //the light vector is directed toward the view point (which is the source of the light)
        // another light vector would be possible
        VectorMath.setVector(lightVector,rayVector[0],rayVector[1],rayVector[2]);

        //Initialization of the colors as floating point values
        double r,g,b;
        r=g=b=0.0;
        double alpha=0.0;
        double opacity=0;

        // Define the increment in each step
        double[]increments=new double[3];
        VectorMath.setVector(increments,rayVector[0]*sampleStep,rayVector[1]*sampleStep,rayVector[2]*sampleStep);

        // Define the distance between the first point and the last point
        double distance=VectorMath.distance(entryPoint,exitPoint);
        int nrSamples=1+(int)Math.floor(distance/sampleStep);

        // Initialize with the first point as the start point
        double[]currentPos=new double[3];
        VectorMath.setVector(currentPos,entryPoint[0],entryPoint[1],entryPoint[2]);

        // Initialization of the colors as floating point values
        TFColor voxel_color=new TFColor(0,0,0,1);
        TFColor colorAux=new TFColor();

        if(compositingMode){
            // 1D transfer function
            do{
                // get the color of the current position by linear interpolation
                colorAux=tFunc.getColor((int)volume.getVoxelLinearInterpolate(currentPos));
                VoxelGradient currentGradient=gradients.getGradient(currentPos);

                // add shading and make it more smooth
                if(shadingMode&&currentGradient.mag>3){
                    colorAux=computePhongShading(colorAux,currentGradient,lightVector,rayVector);
                }

                // compute the voxel by accumulation
                // front-to-back compositing
                voxel_color.r=colorAux.a*colorAux.r+(1-colorAux.a)*voxel_color.r;
                voxel_color.g=colorAux.a*colorAux.g+(1-colorAux.a)*voxel_color.g;
                voxel_color.b=colorAux.a*colorAux.b+(1-colorAux.a)*voxel_color.b;
                voxel_color.a=colorAux.a+(1-colorAux.a)*voxel_color.a;
                for(int i=0;i< 3;i++){
                    currentPos[i]+=increments[i];
                }
                nrSamples--;
            }while(nrSamples>0);

            // save the final result
            r=voxel_color.r;
            g=voxel_color.g;
            b=voxel_color.b;
            alpha=1;

        }

        if(tf2dMode){
            // 2D transfer function
            // get basic intensity and radius from tFunc2D
            short baseIntensity=tFunc2D.baseIntensity;
            double baseRadius=tFunc2D.radius;

            do{
                // get current voxel value and gradient
                double value=volume.getVoxelLinearInterpolate(currentPos);
                VoxelGradient gradinet=gradients.getGradient(currentPos);
                double max=gradinet.mag;

                // get basic color from tFunc2D
                colorAux.r=tFunc2D.color.r;
                colorAux.g=tFunc2D.color.g;
                colorAux.b=tFunc2D.color.b;
                alpha=computeOpacity2DTF(baseIntensity,baseRadius,value,max);
                colorAux.a=alpha;

                // add shading
                if(shadingMode&&gradinet.mag>0){
                    colorAux=computePhongShading(colorAux,gradinet,lightVector,rayVector);
                }

                // compute the color from front to back
                voxel_color.r=alpha*colorAux.r+(1-alpha)*voxel_color.r;
                voxel_color.g=alpha*colorAux.g+(1-alpha)*voxel_color.g;
                voxel_color.b=alpha*colorAux.b+(1-alpha)*voxel_color.b;
                voxel_color.a=alpha+(1-alpha)*voxel_color.a;

                for(int i=0;i< 3;i++){
                    currentPos[i]+=increments[i];
                }
                nrSamples--;
            }while(nrSamples>0);

            // save the result
            r=voxel_color.r;
            g=voxel_color.g;
            b=voxel_color.b;
            alpha=1;

        }

        //computes the color
        int color=computeImageColor(r,g,b,alpha);
        return color;
}


////// VERY IMPORTANT: BISECTION_ACCURACY IS EMBEDDED IN traceRayIso
public int traceRayIso(double[]entryPoint,double[]exitPoint,double[]rayVector,double sampleStep){

        double[]lightVector=new double[3];
        //We define the light vector as directed toward the view point (which is the source of the light)
        // another light vector would be possible
        VectorMath.setVector(lightVector,rayVector[0],rayVector[1],rayVector[2]);

        // To be Implemented

        //Initialization of the colors as floating point values
        double r,g,b;
        r=g=b=0.0;
        double alpha=0.0;
        double opacity=0;


        // To be Implemented this function right now just gives back a constant color
        double[]currentPoint=new double[3];

        //calculate number of sample steps for loop termination
        double steps=VectorMath.distance(entryPoint,exitPoint)/sampleStep;

        double dx=exitPoint[0]-entryPoint[0];
        double dy=exitPoint[1]-entryPoint[1];
        double dz=exitPoint[2]-entryPoint[2];
        double diffsum=Math.sqrt(dx*dx+dy*dy+dz*dz);

        double isovalue=0;

        float voxel=0;float previousvalue=0;

        //loop over all steps
        VoxelGradient gradient=new VoxelGradient();
        for(double step=0;step<steps; step++){
            //save previous point for bisection
            double[]previouspoint=currentPoint.clone();

            //increment current point by 1 step
            currentPoint[0]=entryPoint[0]+step*sampleStep*dx/diffsum;
            currentPoint[1]=entryPoint[1]+step*sampleStep*dy/diffsum;
            currentPoint[2]=entryPoint[2]+step*sampleStep*dz/diffsum;

            //save previous voxel for bisection
            previousvalue=voxel;

            //obtain current voxel value with either linear or cubic interpolation
            //  voxel =  volume.getVoxelLinearInterpolate(currentPoint);
            voxel=volume.getVoxelTriCubicInterpolate(currentPoint);

            //has no access to local variables so we simply replace the function call by the function itself.

            // bisection_accuracy (currentPoint, previouspoint,sampleStep,  previousvalue,voxel,  iso_value, gradient,isovalue);

            //only engage bisection algorithm if we know the iso value is hit somewhere between previous and current voxel
            if(previousvalue<iso_value &&voxel>=iso_value){
                //algorithm tolerance;
                double tol=.01;
                int n=0;
                double[]midpoint=currentPoint.clone();
                //max 10 iterations
                while(n< 10){
                    //midpoint is half the sum of previous and current point
                    midpoint=VectorMath.multiply(VectorMath.add(currentPoint,previouspoint),.5);

                    //get interpolated value at midpoint with either cubic or linear interpolation
                    // double midval = volume.getVoxelLinearInterpolate(midpoint);
                    double midval=volume.getVoxelTriCubicInterpolate(midpoint);

                    //we break the algorithm if the current midpoint interpolation is close enough to the desired iso value
                    if(Math.abs(iso_value-midval)<tol){
                        break;
                    }

                    n++;
                    // here we decide whether we look for the new aproximation to the right or to the left of current midpoint.
                    if(midval<iso_value){
                        previouspoint=midpoint.clone();
                    }else{
                        currentPoint=midpoint.clone();
                    }

                }

                //loop has terminated so midpoint is now our best approximation of isovalue location
                gradient=gradients.getGradient(midpoint);
                isovalue=1;

            }
            //we are done if we have a hit
            if(isovalue>0){
                break;
            }

        }

        //compute phong shading color
        TFColor Phongcolor=computePhongShading(isoColor,gradient,lightVector,rayVector);

        // isoColor contains the isosurface color from the interface
        r=Phongcolor.r;
        g=Phongcolor.g;
        b=Phongcolor.b;


        alpha=isovalue;
        //computes the color
        int color=computeImageColor(r,g,b,alpha);

        return color;
}

// Compute Phong Shading given the voxel color (material color), the gradient, the light vector and view vector
public TFColor computePhongShading(TFColor voxel_color,VoxelGradient gradient,double[]lightVector,
        double[]rayVector){

        //define constants
        double Ka=.1;
        double Kd=.7;
        double Ks=.2;
        int alpha=100;

        //define vectors necessary for Phong Shading. Note that they must be normalized.
        double[]L=VectorMath.normalize(lightVector);
        double[]N={gradient.x,gradient.y,gradient.z};

        N=VectorMath.normalize(N);
        double[]V=VectorMath.normalize(rayVector);
        double factor=2*VectorMath.dotproduct(L,N);
        // reflection vector
        double[]R={factor*N[0]-L[0],factor*N[1]-L[1],factor*N[2]-L[2]};
        R=VectorMath.normalize((R));


        TFColor Lightcolor=new TFColor(1.0,1.0,1.0,1);

        //apply phong model component wise for R,G and B
        double outputR=Lightcolor.r*voxel_color.r*Ka+Lightcolor.r*voxel_color.r*Kd*(VectorMath.dotproduct(N,L))+Lightcolor.r*voxel_color.r*Ks*Math.pow(VectorMath.dotproduct(R,V),alpha);
        double outputG=Lightcolor.g*voxel_color.g*Ka+Lightcolor.g*voxel_color.g*Kd*(VectorMath.dotproduct(N,L))+Lightcolor.g*voxel_color.g*Ks*Math.pow(VectorMath.dotproduct(R,V),alpha);
        double outputB=Lightcolor.b*voxel_color.b*Ka+Lightcolor.b*voxel_color.b*Kd*(VectorMath.dotproduct(N,L))+Lightcolor.b*voxel_color.b*Ks*Math.pow(VectorMath.dotproduct(R,V),alpha);

        TFColor color=new TFColor(outputR,outputG,outputB,1);

        return color;
}


///////////// In GradientVolume
private float interpolate(float g0,float g1,float factor){
        float result=(1-factor)*g0+factor*g1;
        return result;
        }
public  VoxelGradient getGradient(double[]coord){

        if(coord[0]< 0||coord[0]>(dimX-2)||coord[1]< 0||coord[1]>(dimY-2)
        ||coord[2]< 0||coord[2]>(dimZ-2)){
        return new VoxelGradient(0,0,0);
        }
        /* notice that in this framework we assume that the distance between neighbouring voxels is 1 in all directions*/
        int x=(int)Math.floor(coord[0]);
        int y=(int)Math.floor(coord[1]);
        int z=(int)Math.floor(coord[2]);

        float fac_x=(float)coord[0]-x;
        float fac_y=(float)coord[1]-y;
        float fac_z=(float)coord[2]-z;

        VoxelGradient interpolated=new VoxelGradient();

        float t0=interpolate(getGradient(x,y,z).x,getGradient(x+1,y,z).x,fac_x);
        float t1=interpolate(getGradient(x,y+1,z).x,getGradient(x+1,y+1,z).x,fac_x);
        float t2=interpolate(getGradient(x,y,z+1).x,getGradient(x+1,y,z+1).x,fac_x);
        float t3=interpolate(getGradient(x,y+1,z+1).x,getGradient(x+1,y+1,z+1).x,fac_x);
        float t4=interpolate(t0,t1,fac_y);
        float t5=interpolate(t2,t3,fac_y);
        float t6=interpolate(t4,t5,fac_z);


        interpolated.x=t6;

        t0=interpolate(getGradient(x,y,z).y,getGradient(x+1,y,z).y,fac_x);
        t1=interpolate(getGradient(x,y+1,z).y,getGradient(x+1,y+1,z).y,fac_x);
        t2=interpolate(getGradient(x,y,z+1).y,getGradient(x+1,y,z+1).y,fac_x);
        t3=interpolate(getGradient(x,y+1,z+1).y,getGradient(x+1,y+1,z+1).y,fac_x);
        t4=interpolate(t0,t1,fac_y);
        t5=interpolate(t2,t3,fac_y);
        t6=interpolate(t4,t5,fac_z);


        interpolated.y=t6;

        t0=interpolate(getGradient(x,y,z).z,getGradient(x+1,y,z).z,fac_x);
        t1=interpolate(getGradient(x,y+1,z).z,getGradient(x+1,y+1,z).z,fac_x);
        t2=interpolate(getGradient(x,y,z+1).z,getGradient(x+1,y,z+1).z,fac_x);
        t3=interpolate(getGradient(x,y+1,z+1).z,getGradient(x+1,y+1,z+1).z,fac_x);
        t4=interpolate(t0,t1,fac_y);
        t5=interpolate(t2,t3,fac_y);
        t6=interpolate(t4,t5,fac_z);


        interpolated.z=t6;

        return interpolated;

}


public void raycast(double[]viewMatrix){

        //data allocation
        double[]viewVec=new double[3];
        double[]uVec=new double[3];
        double[]vVec=new double[3];
        double[]pixelCoord=new double[3];
        double[]entryPoint=new double[3];
        double[]exitPoint=new double[3];

        // increment in the pixel domain in pixel units
        int increment=1;
        // sample step in voxel units
        int sampleStep=1;
        // reset the image to black
        resetImage();

        // vector uVec and vVec define the view plane,
        // perpendicular to the view vector viewVec which is going from the view point towards the object
        // uVec contains the up vector of the camera in world coordinates (image vertical)
        // vVec contains the horizontal vector in world coordinates (image horizontal)
        getViewPlaneVectors(viewMatrix,viewVec,uVec,vVec);


        // The result of the visualization is saved in an image(texture)
        // we update the vector according to the resolution factor
        // If the resolution is 0.25 we will sample 4 times more points.
        for(int k=0;k<3;k++)
        {
            uVec[k]=res_factor*uVec[k];
            vVec[k]=res_factor*vVec[k];
        }

        // We get the size of the image/texture we will be puting the result of the
        // volume rendering operation.
        int imageW=image.getWidth();
        int imageH=image.getHeight();

        int[]imageCenter=new int[2];
        // Center of the image/texture
        imageCenter[0]=imageW/2;
        imageCenter[1]=imageH/2;

        // imageW/ image H contains the real width of the image we will use given the resolution.
        //The resolution is generated once based on the maximum resolution.
        imageW=(int)(imageW*((max_res_factor/res_factor)));
        imageH=(int)(imageH*((max_res_factor/res_factor)));

        //The rayVector is pointing towards the scene
        double[]rayVector=new double[3];
        rayVector[0]=-viewVec[0];rayVector[1]=-viewVec[1];rayVector[2]=-viewVec[2];

        // compute the volume center
        double[]volumeCenter=new double[3];
        computeVolumeCenter(volumeCenter);

        //compute loop size
        int size=((imageCenter[1]+imageH/2)-(imageCenter[1]-imageH/2))/increment;

        //initialize array of threads equal to size, 1 thread = 1 column of pixels
        threadtest[]threadarray=new threadtest[size];

        int k=0;
        // ray computation for each pixel
        for(int j=imageCenter[1]-imageH/2;j<imageCenter[1]+imageH/2;j+=increment){
            //define thread input parameters
            threadobject params=new threadobject(viewVec,uVec,vVec,increment,sampleStep,imageW,imageCenter[0],rayVector,j);
            // initialize thread
            threadarray[k]=new threadtest(params);
            //tell thread to start executing
            threadarray[k].start();


            k++;
        }
        //after initializing all threads, main thread must wait for all threads to terminate
        for(threadtest thread:threadarray){
            try{
                thread.join();
            }
                catch(InterruptedException e){}
            }
        }

public class threadtest extends Thread {
    private threadobject params;

    private threadtest(threadobject params) {
        this.params = params;
    }

    @Override
    public void run() {
        int i1 = params.getI1();
        int imageW = params.getImageW();
        int increment = params.getIncrement();
        double[] viewVec = params.getViewVec();
        double[] uVec = params.getuVec();
        double[] vVec = params.getvVec();
        double[] rayVector = params.getRayVector();
        int sampleStep = params.getSampleStep();
        int j = params.getJ();


        //this is basically a copy of the original code
        for (int i = i1 - imageW / 2; i < i1 + imageW / 2; i += increment) {
            double[] pixelCoord = new double[3];
            double[] entryPoint = new double[3];
            double[] exitPoint = new double[3];
            // compute starting points of rays in a plane shifted backwards to a position behind the data set
            computePixelCoordinatesBehindFloat(pixelCoord, viewVec, uVec, vVec, i, j);
            // compute the entry and exit point of the ray
            computeEntryAndExit(pixelCoord, rayVector, entryPoint, exitPoint);
            if ((entryPoint[0] > -1.0) && (exitPoint[0] > -1.0)) {
                int val = 0;
                if (compositingMode || tf2dMode) {
                    val = traceRayComposite(entryPoint, exitPoint, rayVector, sampleStep);
                } else if (mipMode) {
                    val = traceRayMIP(entryPoint, exitPoint, rayVector, sampleStep);
                } else if (isoMode) {
                    val = traceRayIso(entryPoint, exitPoint, rayVector, sampleStep);
                }
                for (int ii = i; ii < i + increment; ii++) {
                    for (int jj = j; jj < j + increment; jj++) {
                        image.setRGB(ii, jj, val);
                    }
                }
            }

        }
    }
}

package tudelft.cgv.volvis;

public class threadobject {

    //data class to pass parameters to threads
    private final double[] viewVec;
    private final double[] uVec;
    private final double[] vVec;
    private final int increment;
    private final int sampleStep;
    private final int imageW;
    private final int i1;
    private final double[] rayVector;
    private final int j;

    public threadobject(double[] viewVec, double[] uVec, double[] vVec, int increment, int sampleStep, int imageW, int i1, double[] rayVector, int j) {
        this.viewVec = viewVec;
        this.uVec = uVec;
        this.vVec = vVec;
        this.increment = increment;
        this.sampleStep = sampleStep;
        this.imageW = imageW;
        this.i1 = i1;
        this.rayVector = rayVector;
        this.j = j;
    }

    ///////////////////  VectorMath.java
    public static double[] normalize(double[] v) {
        double[] result = new double[3];
        // uses vectormath.length, just to be clear it is not the length of the array (3)
        if (length(v) > 0) {
            result[0] = v[0] / length(v);
            result[1] = v[1] / length(v);
            result[2] = v[2] / length(v);

            return result;
        } else {
            return v;
        }
    }

    public static double[] add(double[] v, double[] w) {

        double[] result = new double[v.length];
        for (int i = 0; i < v.length; i++) {
            result[i] = v[i] + w[i];
        }
        return result;
    }

    public static double[] multiply(double[] v, double w) {

        double[] result = new double[v.length];
        for (int i = 0; i < v.length; i++) {
            result[i] = v[i] * w;
        }
        return result;
    }

    public double[] getViewVec() {
        return viewVec;
    }

    public double[] getuVec() {
        return uVec;
    }

    public double[] getvVec() {
        return vVec;
    }

    public int getIncrement() {
        return increment;
    }

    public int getSampleStep() {
        return sampleStep;
    }

    public int getImageW() {
        return imageW;
    }
}

    public int getI1() {
        return i1;
    }

    public double[] getRayVector() {
        return rayVector;
    }

    public int getJ() {
        return j;
    }