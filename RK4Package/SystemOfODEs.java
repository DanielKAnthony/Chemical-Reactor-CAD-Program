package RK4Package;

public class SystemOfODEs extends reactors
{
  /*this class has an array of pfr ODEs that make up the system of ODEs*/
  
  private PFR[] odeSystem;
  
  public SystemOfODEs(PFR[] odeSystem)
  {
    super();
    if(odeSystem==null)
      System.exit(0);
    for(int i=0;i<odeSystem.length;i++)
    {
      if(odeSystem[i]==null)
        System.exit(0);
    }
    
    this.odeSystem = new PFR[odeSystem.length];
    for(int i=0;i<odeSystem.length;i++)
    {
      this.odeSystem[i]=odeSystem[i].clone();
      
    }
  }//end of constructor
  
  public SystemOfODEs(SystemOfODEs source)
  {
    super(source);
    if(source==null)
      System.exit(0);
    
    this.odeSystem = new PFR[source.odeSystem.length];
    for(int i=0;i<source.odeSystem.length;i++)
    {
      this.odeSystem[i]=source.odeSystem[i].clone();
    }
  }//end of copy constructor
  
  public SystemOfODEs clone()
  {
    return new SystemOfODEs(this);
  }//end of clone
  
  public PFR[] getOdeSystem()
  {
    PFR[] copy =new PFR[this.odeSystem.length];
    for(int i=0;i<this.odeSystem.length;i++)
      copy[i]=this.odeSystem[i].clone();
    return copy;
  }//end of accessor
  
  public boolean setOdeSystem(PFR[] odeSystem)
  {
    if(odeSystem==null || odeSystem.length!=this.odeSystem.length)
      return false;
    for(int i=0;i<odeSystem.length;i++)
    {
      if(odeSystem[i]==null)
        return false;
    }
    this.odeSystem = new PFR[odeSystem.length];
    for(int i=0;i<odeSystem.length;i++)
      this.odeSystem[i]=odeSystem[i].clone();
    return true;
  }//end of mutator 
  
  public double[] solveFlowRates(double[] F, double V)
  {
    /*this method calculates the flow rates of all species at the outlet of the reactor volume by calling 
     * the RK4Fel ODE solver method and passing it the system of ODEs*/
    double delx=V/1000;//the step size is set to the volume size divided by 1000
    double[] flow=new double[F.length];
    flow= RK4Fel.solveODE(0., F, V, delx, 10001, this.odeSystem);
    this.reset();
    
    return flow;
  }
    
  }
//end of class