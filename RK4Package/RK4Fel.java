package RK4Package;

public class RK4Fel
{    
  
  public static double[] solveODE(double x_0, double[] y_0, double x_f, double delx, int maxIt, NumMethod[] f) 
  {
    double x=x_0;//initialize x
    double[] y= new double[y_0.length];
    for(int i=0;i<y_0.length;i++)
      y[i] = y_0[i];//initialize y
    int i=0;
    double[] error = new double[y_0.length];
    double tol = 0.001;//tolerance is set to 0.001
    double[] y_next = new double[y_0.length];
    double[] z_next = new double[y_0.length];
    double[] s = new double[y_0.length];
    boolean check = true;
    double s_min = 10.;
    int counter=0;
    
    while(x<=x_f && i<maxIt)
    {
      do
      {
      
      double[] k1 = new double[y_0.length];
      double[] k2 = new double[y_0.length];
      double[] k3 = new double[y_0.length];
      double[] k4 = new double[y_0.length];
      double[] k5 = new double[y_0.length];
      double[] k6 = new double[y_0.length];
      //calculating the k1 values
      for(int j=0;j<y_0.length;j++)
      {
        k1[j]=delx*f[j].calculateVals(x, y);
      }
      //calculating the k2 values
      for(int j=0;j<y_0.length;j++)
      {
        double[] yTemp = new double[y.length];
        for(int l=0; l<y.length;l++)
          yTemp[l]=y[l]+k1[l]/4.;
        k2[j]=delx*f[j].calculateVals(x+delx/4.,yTemp);
      }
      //calculating the k3 values
      for(int j=0;j<y_0.length;j++)
      {
        double[] yTemp = new double[y.length];
        for(int l=0; l<y.length;l++)
          yTemp[l]=y[l]+3.*k1[l]/32.+9.*k2[l]/32.;
        k3[j]=delx*f[j].calculateVals(x+3.*delx/8.,yTemp);
      }
      //calculating the k4 values
      for(int j=0;j<y_0.length;j++)
      {
        double[] yTemp = new double[y.length];
        for(int l=0; l<y.length;l++)
          yTemp[l]=y[l]+1932.*k1[l]/2197.-7200.*k2[l]/2197.+7296.*k3[l]/2197.;
        k4[j]=delx*f[j].calculateVals(x+12.*delx/13.,yTemp);
      }
      //calculating the k5 values
      for(int j=0;j<y_0.length;j++)
      {
        double[] yTemp = new double[y.length];
        for(int l=0; l<y.length;l++)
          yTemp[l]=y[l]+439.*k1[l]/216.-8*k2[l]+3680.*k3[l]/513.-845.*k4[l]/4104.;
        k5[j]=delx*f[j].calculateVals(x+delx,yTemp);
      }
      //calculating the k6 values
      for(int j=0;j<y_0.length;j++)
      {
        double[] yTemp = new double[y.length];
        for(int l=0; l<y.length;l++)
          yTemp[l]=y[l]-8.*k1[l]/27.+2.*k2[l]-3544.*k3[l]/2565.+1859.*k4[l]/4104.-11.*k5[l]/40.;
        k6[j]=delx*f[j].calculateVals(x+delx/2.,yTemp);
      }
      //calculating next y values
      for(int j=0;j<y_0.length;j++)
      {
        y_next[j]=y[j]+25.*k1[j]/216.+1408.*k3[j]/2565.+2197.*k4[j]/4104.-k5[j]/5.;
      }
      //calculating next z values
      for(int j=0;j<y_0.length;j++)
      {
        z_next[j]=y[j]+16.*k1[j]/135.+6656.*k3[j]/12825.+28561.*k4[j]/56430.-9.*k5[j]/50.+2.*k6[j]/55.;
      }
      //calculating the error
      for(int j=0;j<y_0.length;j++)
      {
        error[j] = (1./delx)*Math.abs(z_next[j]-y_next[j]);
      }
      //calculating the s parameter
      for(int j=0;j<y_0.length;j++)
      {
        s[j] = 0.84*Math.pow(tol/error[j],0.25);
      }
      //checking if the error is less than the tolerance
      for(int j=0;j<y_0.length;j++)
      {
        if(error[j]>=tol)
          check = false;
      }
      //calculate the minimum s parameter
      for(int j=0;j<y_0.length;j++)
      {
        if(s[j]<=s_min)
          s_min = s[j];
      }
      //changing the step size
      delx = s_min*delx;
      counter++;
      }while(check==false && counter<=10000);

      //updating x and y values
      x=x+delx;
      for(int j=0;j<y_0.length;j++)
      {
        y[j] = y_next[j];
      }
      i++;
      
    }    
    if(i>maxIt) System.out.println("endpoint not reached within maximum specified initerations in RK4 integrate method");
    
    return y;//returning the exiting flow rates of all species at the end of the reactor volume
  }//end of solveODE method
}//end of RK4 class
