package RK4Package;

//the instance variable g_r will be set equal to the reaction rate calculated in RateInTermsOfFlow

public abstract class reactors {
  
  private rxn g_r;
  protected double[][] g_stoich;
  protected double[] g_k;
  protected double[] g_keq;
  protected int[] g_base;
  protected double[] g_F;
  protected double[] g_kd;
  protected double g_ct0;
  protected int g_numofrxn;
  protected double[][] g_rxnorder;
  protected boolean g_phase;
  protected double[] g_initflows;
  public reactors() {
   
    //set global variables to null of zero initially
   this.g_r=null;
   this.g_stoich=null;
   this.g_k=null;
   this.g_keq=null;
   this.g_base=null;
   this.g_F=null;
   this.g_kd=null;
   this.g_ct0=0;
   this.g_numofrxn=0;
   this.g_rxnorder=null;
   this.g_phase=false;
   this.g_initflows=null;
       
  }
  
  //set global variables to null or zero
  public reactors(reactors source) {
   if(source==null) System.exit(0);
   this.g_r=source.g_r.clone(); 
   this.g_stoich=null;
   this.g_k=null;
   this.g_keq=null;
   this.g_base=null;
   this.g_F=null;
   this.g_kd=null;
   this.g_ct0=0;
   this.g_numofrxn=0;
   this.g_rxnorder=null;
      this.g_phase=false;
   this.g_initflows=null;
  }
  
  //resets globall variables to null or zero
  protected void reset() {
   
    this.g_r=null; 
   this.g_stoich=null;
   this.g_k=null;
   this.g_keq=null;
   this.g_base=null;
   this.g_F=null;
   this.g_kd=null;
   this.g_ct0=0;
   this.g_numofrxn=0;
   this.g_rxnorder=null;
      this.g_phase=false;
   this.g_initflows=null;
  
  }

public abstract reactors clone(); //now available for each child class

protected boolean setR(rxn g_r) {
  if(g_r==null) return false;
  this.g_r=g_r.clone();
  return true;
}

//methods available to be called by PFR

protected double returnReactionRate(double []F, double [][]stoich, double[] k, int numofrxn, double ct0, double[] keq, int[] base, int p, double[][]rxnorder, boolean phase, double[] initflows) {
 if(this.g_r==null) System.exit(0);
  return this.g_r.calcRateInTermsFlow(F, stoich, k, numofrxn, ct0, keq, base, p, rxnorder, phase, initflows); 
}

protected double returnReactionRate(double []F, double [][]stoich, double[] k, int numofrxn, double ct0, int[] base, int p, double[][]rxnorder, boolean phase, double[] initflows) {
 if(this.g_r==null) System.exit(0);
  return this.g_r.calcRateInTermsFlow(F, stoich, k, numofrxn, ct0, base, p, rxnorder, phase, initflows); 
}
  
}