package RK4Package;

public class PFR extends reactors implements NumMethod, Values{
 
 private int cntr;//represents the array index for system of ODEs
 
 public PFR(int cntr)
 {
   super();
   this.cntr=cntr;
 }//end of constructor
 
 public PFR(PFR source) {
   
  this.cntr=source.cntr; 
 }
  
  public PFR clone()
 {
   return this;
 }//end of clone
  //sets all global variables to values inputted by user or retrieved from text files
  public void setVals(double[] F, double[][] stoich, double [] k, double [] keq, int[] base, double ct0, int numofrxn, rxn r, double V, double[] kd, double[][]rxnorder, boolean phase, double[] initflows) {
    super.setR(r);
    
    super.g_stoich= new double[stoich.length][stoich[0].length];
    for(int i=0; i<stoich.length;i++) {
    for(int j=0; j<stoich[0].length;j++) 
     super.g_stoich[i][j]=stoich[i][j];
    }
    super.g_rxnorder= new double[rxnorder.length][rxnorder[0].length];
    for(int i=0; i<rxnorder.length;i++) {
    for(int j=0; j<rxnorder[0].length;j++) 
     super.g_rxnorder[i][j]=rxnorder[i][j];
    }

    super.g_initflows=new double[initflows.length];
    for(int i=0; i<initflows.length; i++)
    super.g_initflows[i]=initflows[i];
    
    super.g_k=new double[k.length];
    for(int i=0; i<k.length; i++)
    super.g_k[i]=k[i];
    
    super.g_keq=new double[keq.length];
    for(int i=0; i<keq.length; i++)
      super.g_keq[i]=keq[i];
    
    super.g_base=new int[base.length];
    for(int i=0; i<base.length;i++) 
      super.g_base[i]=base[i];
    
    super.g_F=new double[F.length];
    for(int i=0; i<F.length;i++) 
      super.g_F[i]=F[i];
    
    super.g_kd=new double[kd.length];
    for(int i=0; i<kd.length;i++) 
      super.g_kd[i]=kd[i];
    
    super.g_ct0=ct0;
    super.g_numofrxn=numofrxn;
    super.g_phase=phase;
  }
  //calculates the diffusion through the membrane
    public double[] calcTransportR(double[] F) {
    double sum=0;
     for(int i=0; i<F.length; i++) //calculates sum of flow rates
    {
      sum+=F[i];
    }
  
    double[] R=new double[F.length];
    for(int i=0; i<F.length; i++)
    R[i]=g_kd[i]*g_ct0*F[i]/sum;//diffusion equals the diffusion constant times the concentration of the species diffusing through the membrane
  
    return R;
  }
  
//this method is called by the RK4Fel solver to obtain the right hand side value of the ODE
  public double calculateVals(double x, double[] y) {
    //x represents V, y represents flow
    
    double[] rtot=new double[y.length];
    double[] r=new double[g_numofrxn];
    double[][] rall=new double[g_numofrxn][y.length];
    double[] R=new double[y.length];
    
    for(int p=0; p<g_numofrxn;p++) {
      if(g_keq[p]==0)
    r[p]=returnReactionRate(y, g_stoich, g_k, g_numofrxn, g_ct0, g_base, p, g_rxnorder, g_phase, g_initflows);
      if(g_keq[p]!=0)
    r[p]=returnReactionRate(y, g_stoich, g_k, g_numofrxn, g_ct0, g_keq, g_base, p, g_rxnorder, g_phase, g_initflows);
    }//calls returnReactionRate from reactors
    
    R=calcTransportR(y);//calls transport rate method
    
    for(int p=0; p<g_numofrxn;p++) {
      for(int i=0; i<g_F.length;i++) {
       rall[p][i]=r[p]*g_stoich[p][i]/g_stoich[p][g_base[p]]; //p=rxn, i=components
      }
    }

     for(int i=0; i<g_F.length;i++) { 
      rtot[i]=0;
      for(int p=0; p<g_numofrxn; p++) { 
    rtot[i]+=rall[p][i]; 
    
      }
      rtot[i]=rtot[i]-R[i];
     }//calculates the total reaction rate for each species (right hand side of ODE)
    
return rtot[this.cntr];
    
    
  }
}
  
  
  