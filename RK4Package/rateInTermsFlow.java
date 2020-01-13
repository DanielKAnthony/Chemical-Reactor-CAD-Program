package RK4Package;

public class rateInTermsFlow extends rxn {

  public rateInTermsFlow clone() {
    return this;
  }
  

  //calculates reaction rate for reversible reactions
  public double calcRateInTermsFlow(double []F, double [][]stoich, double[] k, int numofrxn, double ct0, double[] keq, int[] base, int p, double[][]rxnorder, boolean phase, double[] initflows) {

    double reactant=1;//product of reactant concentration terms
    double product=1;//product of product concentration terms
    double sum=0;
    double r=0;
    double[][] rall=new double[numofrxn][F.length];//reaction rates for each species per reaction
    double[] rtot=new double[F.length];//sum of all reaction rates corresponding to each species

    if(phase){
      for(int a=0; a<initflows.length; a++) {
        sum+=initflows[a]; 
      }
    }//the total flow for liquid phase reactions is constant and equal to the initial total flow
      
      else{ for(int i=0; i<F.length; i++) //calculates sum of flow rates for gas only
    {
      sum+=F[i];
      }}//the total flow for gas phase reactions is changing and must be adjusted each time 
    
      reactant=1;
      product=1;
  for (int i=0; i<F.length; i++) //for loop for each components
  { 
    
    if(stoich[p][i]<0)
  reactant *= Math.pow(ct0*F[i]/sum,Math.abs(rxnorder[p][i]));
   else if(stoich[p][i]>0) product*= Math.pow(ct0*F[i]/sum, (rxnorder[p][i]));
   //calculates the product of reactant concentration terms
}

  r= -k[p]*reactant+(k[p]/keq[p])*product;
  //calculates the reaction rate for base species
  
    
    
    return r;
  }

  //overloads previous method to calculate rate for irreversible reactions
  public double calcRateInTermsFlow(double []F, double [][]stoich, double[] k, int numofrxn, double ct0, int[] base, int p, double[][]rxnorder, boolean phase, double[] initflows) {
    double reactant=1;
    double sum=0;
    double r=0;
    double[][] rall=new double[numofrxn][F.length];
    double[] rtot=new double[F.length];
    if(phase){
          for(int a=0; a<initflows.length; a++) {
        sum+=initflows[a]; 
      }
    }
      
      else{
    
    for(int i=0; i<F.length; i++) //calculates sum of flow rates
    {
      sum+=F[i];
    }
      }

      reactant=1;
  for (int i=0; i<F.length; i++) //for loop for each components
  { 
    
    if(stoich[p][i]<0)
  reactant *= Math.pow(ct0*F[i]/sum,Math.abs(rxnorder[p][i]));
   
}

  r= -k[p]*reactant;
  
  //formulas for irreversible reactions
  
    return r;
    
      }
      
} 
    



  
  