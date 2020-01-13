package RK4Package;

public abstract class rxn {
  /*This abstract class set out three methods that every child object of this class such as rateInTermFlow must create*/
  
  public abstract double calcRateInTermsFlow(double []F, double [][]stoich, double[] k, int numofrxn, double ct0, double[] keq, int[] base, int p, double[][]rxnorder, boolean phase, double[] initflows);
  public abstract double calcRateInTermsFlow(double []F, double [][]stoich, double[] k, int numofrxn, double ct0, int[] base, int p, double[][]rxnorder, boolean phase, double[] initflows);  
  public abstract rxn clone();
}