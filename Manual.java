import RK4Package.*;
import PFRPackage.*;
import java.util.*;
import java.lang.Math;

public class Manual{

 private static double safeGetDouble(String message) {
  while (true) {
   try {
    if (!message.equals(""))
     System.out.println(message);
    Scanner scan = new Scanner(System.in);
    return scan.nextDouble();
   } catch (Exception e) {
    System.out.println("That isn't a number!");
   }
  }
    
 }

 private static int safeGetInt(String message) {
  while (true) {
   try {
    if (!message.equals(""))
     System.out.println(message);
    Scanner scan = new Scanner(System.in);
    return scan.nextInt();
   } catch (Exception e) {
    System.out.println("An integer is needed");
   }
  }
    
 }

 private static char safeGetChar(String message) {
  while (true) {
   try {
    if (!message.equals(""))
     System.out.println(message);
    Scanner scan = new Scanner(System.in);
 char c = scan.next().toUpperCase().charAt(0);
  if ('A' <= c && c <= 'Z') 
   return c;
  else
   System.out.println("Letters only!");
 } catch (Exception e) {
      System.out.println("That isn't valid");
    }
  }
 }

 private static boolean safeGetYesNo (String message) {
 char yesno = 'z';
 while (yesno != 'Y' && yesno != 'N') {
 yesno = safeGetChar(message + " [enter y or n]");
 }
 return yesno == 'Y';
 }

public static void waitForUser(){
  Scanner prompt = new Scanner(System.in);
  System.out.println("\n\nPress ENTER to continue");
  prompt.nextLine();
}

public static void pfrPrompts(){
    Scanner scan = new Scanner(System.in);
    boolean isLiq;
    int pfrnum = safeGetInt("Number of PFRs: ");
    double vol = safeGetDouble("Volume of one PFR: ");
    int rxNum = safeGetInt("Number of Reactions:");
    int specs = safeGetInt("Number of species in the system: ");
    boolean phase_check = true;
    int liq_int = -1;
    while(phase_check){
      liq_int = safeGetInt("Enter 1 for liquid phase or 0 for gas phase: ");
      if(liq_int == 1 || liq_int == 0) phase_check = false;
    }
    
    if(liq_int == 1) isLiq = true;
    else isLiq = false;

    for(int i = 0;i<25;i++){
      System.out.print("_");
    }
    System.out.println("\n\nNOTE: For each reaction, if a species is not present\nenter a coefficient of 0. If a species is a product, enter its\ncoefficient as a negative value.\n");
    for(int i = 0;i<25;i++){
      System.out.print("_");
    }

    waitForUser();

    char specNames = 'A';
    double[][] spec_coeffs = new double[rxNum][specs];
    double[] flows = new double[specs];

    for(int i = 0;i<rxNum;i++){
      System.out.println("For reaction "+(i+1)+":");
      for(int j = 0;j<specs;j++){
        spec_coeffs[i][j] = safeGetDouble("Coefficient of "+specNames+":");
        specNames++;
      }
      specNames = 'A';
    }

    System.out.println("Next, please enter inital flow rates.\nEnter 0 if none.");

    waitForUser();

    for(int i = 0;i<specs;i++){
      flows[i] = safeGetDouble("Initial flow rate of "+specNames+":");
      specNames++;
    }

    System.out.println("\nNext, please enter a k value for each reaction.");
    System.out.println("When prompted to indicate which species the k value is in terms of,\nEnter an integer such that 1 corresponds to A, 2 corresponds to B, etc.");

    waitForUser();

    double[] userK = new double[rxNum];
    int[] k_index = new int[rxNum]; 

    for(int i = 0;i<rxNum;i++){
      System.out.println("For reaction "+(i+1)+":");
      userK[i] = safeGetDouble("Enter reaction constant:");
      k_index[i] = safeGetInt("Which integer corresponds to the species which k is in terms of?");
    }


    double[] kInTermsOf = new double[k_index.length];
    for(int i = 0;i<k_index.length;i++){
      k_index[i] -= 1;
      kInTermsOf[i] = spec_coeffs[i][k_index[i]];
    }

    double[][] rxnK = new double[rxNum][specs];

    for(int i = 0;i<rxNum;i++){
      for(int j = 0;j<specs;j++){
        if((k_index[i]) == j){
          rxnK[i][j] = userK[i];
        }
        else{
          rxnK[i][j] = 0;
        }
      }
    }

    RK4Flows rkflow = new RK4Flows();
    FlowSolver flow_solve = new FlowSolver();

    double[] interstage = new double[] {flows[1],1,1,1,1};
    double[] test = rkflow.iterFlows(isLiq,flows,spec_coeffs,rxNum,rxnK,kInTermsOf,interstage,pfrnum,vol);

    specNames = 'A';
    System.out.println("Exit flows are as follows:\n");
    for(int i = 0;i<test.length;i++){
      System.out.println("F_"+specNames+" = "+test[i]);
      specNames++;
    }
 }

 public static void main()
 {

 boolean flag = true; 
 int multirxn = -1; 
 int numrxn = -1;
 int rxnord= -1;
 double rxnk=-1;
 int rnxk_allotment = -1;
 boolean rxnordercustom = false;
 boolean rxnmembrane = false;
 //nous avons besoin de fichiers d'entree pour dire la phase liquide ou gazeuse et la section io pour gere cette
 boolean phase = safeGetYesNo("Press 'y' for liquid phase, or press 'n' for gas phase.");

    
 Scanner pfrQuestion = new Scanner(System.in);
 String pfr_char = "placeholder";
 boolean pfrFlag = true;
 while(pfrFlag){
  try{
    System.out.println("Is the system a series of PFRs? (y/n)");
    pfr_char = pfrQuestion.nextLine();
    if(pfr_char.charAt(0) == 'y' || pfr_char.charAt(0) == 'n') pfrFlag = false;
    else System.out.println("Please enter y or n");
  }catch(Exception e){
    System.out.println("Invalid");
    pfrQuestion.nextLine();
  }
}

if(pfr_char.charAt(0) == 'y') pfrPrompts();
else{
 
 multirxn = safeGetInt("How many chemical species are present (both reactants and products)?");
  if (multirxn <= 1) {
  System.out.println("one or less species, no computations needed. Exiting.");
  return;
  }

  if(multirxn>27){

  System.out.println("The number of species input is a bit much, try again"); //Because 
  multirxn = safeGetInt("How many chemical species are present (both reactants and products)?");
  }
 numrxn = safeGetInt("How many reaction sets are present?");
 if (numrxn <= 0) {
  System.out.println("zero or less rxns, no computations needed. Exiting.");
  return;
 }


 //SPECIES ARRAY BEING CREATED
 char[] reactname = new char[multirxn]; 
 double[][] reactstoic = new double[numrxn][multirxn];
 char start = 'A';

 System.out.print("Species names are as follows: ");
 for (int i=0; i<reactname.length; i++) {
  System.out.printf("%2c ", start);
  reactname[i]=start;
  start++;
 }
 System.out.println();

 //STOICHIOMETRY ARRAY BEING CREATED
 boolean onepositive;  // need at least one of each
 boolean onenegative; // need at least one of each
 for(int i=0; i<reactstoic.length; i++) {
  do {
  onepositive  = false;
  onenegative = false;
  for (int j=0; j<reactstoic[i].length; j++) {
  reactstoic[i][j]=safeGetDouble("input the species " + reactname[j] +
     " stoichiometric coefficient for reaction set number " +(i+1)+
     " \n(NOTE positive numbers are for products, negatives for reactants)");
  if (reactstoic[i][j] < 0) {
  onenegative = true;
  } else if (reactstoic[i][j] > 0) {
  onepositive = true;
  }
  }
  if (!onepositive || !onenegative) 
  System.out.println("You can't have a reaction without some reactants (negative coefficients) and some products (positive coefficients)");
  } while (! (onepositive && onenegative));
 }


 //GETTING RXN ORDERS
 rxnordercustom = safeGetYesNo("Do you want a custom reaction rate exponents?");

 double[][] rxnorder= new double[numrxn][multirxn];
 if(rxnordercustom == true) {
  for(int i=0; i<reactstoic.length; i++){
  for(int j=0; j<reactstoic[i].length; j++){
  rxnorder[i][j] = safeGetDouble("input the rxn order exponent you want in rxn set " +(i+1)+" for the concentration of species " +reactname[j]);
  }
  }
 } else {    
  rxnord = -1;
  while (rxnord < 0) {
  rxnord = safeGetInt("Press 0 if you want elementary reaction kinetics, else type the order (typing 1 will be first order, 2 will be second order etc.)");
  if (rxnord < 0)
  System.out.println("You cannot have a negative order.");
  }

  if (rxnord == 0){
  for (int i=0; i<reactstoic.length; i++){
  for(int j=0; j<reactstoic[i].length; j++)
  rxnorder[i][j]=reactstoic[i][j];
  }
  }
  else {
  for(int i=0; i<reactstoic.length; i++){
  for(int j=0; j<reactstoic[i].length; j++)
  rxnorder[i][j]=rxnord;
  }
  }
 }

  //GETTING THE Keq
 double[] keq = new double[reactstoic.length];
 for (int i=0; i<keq.length; i++) {
  keq[i]=safeGetDouble("Input the equilibrium constant Keq for reaction set "+(i+1)+ "\nor type 0 if irreversible");
 }


 //GETTING THE kforwards!
 double[] forwardk= new double[numrxn]; //just 1D array
 int[] base = new int[numrxn];
 for (int i=0; i<forwardk.length; i++)
 {
  forwardk[i]=safeGetDouble("Enter the known reaction rate constant for reaction "+(i+1)+" below");
  while (true) {
  int ex = (int) safeGetChar("Enter which species the above constant was for \n(If the rate constant entered above was for species A press A, if it was for species B press B, etc. )");
  ex = ex - 'A';
  if (reactstoic[i][ex] >=0 && keq[i]==0) {
  System.out.println("you can't have a reaction constant in terms of a product in an irreversible reaction");
  System.out.println("for reference,the stoichiometry for your input reactions are:");
  System.out.print("\n    ");
  for(int x=0; x<reactname.length; x++){
  System.out.printf("%8c ", reactname[x]);
  }
    for(int j=0; j<reactstoic.length; j++){
  System.out.printf("\n rxn%-3d",j+1);
  for(int w=0; w<reactstoic[j].length; w++)
   System.out.printf(" %+7.2f ", reactstoic[j][w]);
  System.out.print("\n");}

  continue;
  }

  if (ex < 0 || ex >= reactname.length) {
  System.out.println("You input a species that doesn't exist in this reaction set");
  System.out.println("for reference, the stoichiometry for your input reactions are:");
  System.out.print("\n    ");
  for(int x=0; x<reactname.length; x++){
  System.out.printf("%8c ", reactname[x]);
  }
  for(int j=0; j<reactstoic.length; j++){
  System.out.printf("\n rxn%-3d",j+1);
  for(int w=0; w<reactstoic[j].length; w++)
   System.out.printf(" %+7.2f ", reactstoic[j][w]);
  System.out.print("\n");}

  continue;
  }
  if (reactstoic[i][ex] == 0) {
  System.out.println("you input a species that doesn't exist in this specific reaction");
  System.out.println("for reference,the stoichiometry for your input reactions are:");
  System.out.print("\n    ");
  for(int x=0; x<reactname.length; x++){
  System.out.printf("%8c ", reactname[x]);
  }

  for(int j=0; j<reactstoic.length; j++){
  System.out.printf("\n rxn%-3d",j+1);
  for(int w=0; w<reactstoic[j].length; w++)
   System.out.printf(" %+7.2f ", reactstoic[j][w]);
  System.out.println();
  }
  continue;
  }
  base[i]=ex;
  break;
  }
 }
 
 //Getting parameters for kd for membrane diffusion
 rxnmembrane= safeGetYesNo("If you want a membrane reactor press 'y' now");
 double[] kd = new double[multirxn];
 if(rxnmembrane == true) {
  for (int i=0; i<kd.length; i++) {
  kd[i]=safeGetDouble("Input the diffusivity constant of species "+reactname[i]+" or type 0 if there is no diffusion");
  }
 } else {
  for (int i=0; i<kd.length; i++)
  kd[i]=0;
 }

 //GETTING THE FLOWS
 double[] flow = new double[multirxn];
 for (int i=0; i<flow.length; i++) {
  flow[i] = safeGetDouble("Input inlet flow (in mol/min)to the reactor for species "+reactname[i]);
 }
 double[] initflows = new double[flow.length];
 if(phase==true)    for(int i=0; i<flow.length; i++) //calculates sum of flow rates for liq only
    {
      initflows[i]=flow[i];
    }
 
 //GETTING parameters for the cto
 double v0 = safeGetDouble("please provide the total entering volumetric flow rate in dm^3/min");
   while(v0<=0){
     System.out.println("The volumetric flow you input is invalid, try again");
   v0 = safeGetDouble("please provide the total entering volumetric flow rate in dm^3/min");
   }
   double ct0 = 0;
   for (int i=0; i<flow.length; i++) {
    ct0 = ct0 + flow[i];
   }
   ct0 = ct0/v0;
    
//DECIDING IF WE WANT TO SOLVE FOR VOLUME OR PRODUCTION/SELECTIVITY
System.out.println("If you want to solve with known volume press 'V' , if you want to solve for largest production press 'P', or if you want to solve for highest selectivity press 'S'");
char option; //MAKE SURE OPTION IS ONE OF V P OR S
option = safeGetChar("");
double v=0;

    if(option=='V'){
      
       //GETTING VOLUME 
       v=safeGetDouble("please provide the reactor volume (in dm^3)");                          

       while(v<=0) {
       System.out.println("the input volume is invalid, please try again");
       v=safeGetDouble("please provide the reactor volume (in dm^3)");
       }
         rateInTermsFlow rate= new rateInTermsFlow();
         PFR[] pfr = new PFR[flow.length];

         for(int i=0; i<flow.length; i++) {
           pfr[i]=new PFR(i);
           pfr[i].setVals(flow, reactstoic, forwardk, keq, base, ct0, numrxn, rate, v, kd, rxnorder, phase, initflows);
         }
         SystemOfODEs client = new SystemOfODEs(pfr);


         double[] flow2=client.solveFlowRates(flow, v);
         for(int i=0; i<flow2.length; i++){
           System.out.println( "flow is :" +flow2[i]+"  for species "+reactname[i]);
         } 
     }
                      
                      
   if(option=='P'){
     //iterate through different volumes to get the highest production 
     int tonum=0;
     boolean flagger=true;
     char prodspec='z';
     while(flagger){
     prodspec=safeGetChar("input the species you want the production to be in terms of (it has to be a product)");
     tonum= (int) prodspec-'A';
     for(int z=0; z<reactstoic.length; z++){
       if(reactstoic[z][tonum]>0) flagger=false;}                                         
    }
                             
    double maxFlow;
    double nextFlow;
    double inc = 10000;
    double cap = 100000;
    double loBound = 0;
    double hiBound = cap;
    double maxV = 0;
    int initialMaxV = 0;
    
    do{
      maxFlow = -1;
      for(v = loBound; v <= hiBound; v = v + inc)
      {
        rateInTermsFlow rate= new rateInTermsFlow();
        PFR[] pfr = new PFR[flow.length];
        
        for(int i=0; i<flow.length; i++) {
          pfr[i]=new PFR(i);
          pfr[i].setVals(flow, reactstoic, forwardk, keq, base, ct0, numrxn, rate, v, kd, rxnorder, phase, initflows);
        }
        SystemOfODEs client = new SystemOfODEs(pfr);
        double[] flow2=client.solveFlowRates(flow, v);
        
        nextFlow = flow2[tonum];
        
        if(nextFlow > maxFlow && Math.abs(nextFlow-maxFlow)/inc > 0.001)
        {
          maxFlow = nextFlow;
          maxV = v;
        }
        else
          break;
      }
      if(maxV==0)
        loBound=0;
      else
        loBound = maxV - inc;
      if(maxV==cap)
        hiBound=cap;
      else
        hiBound = maxV + inc;
      inc = inc/2;
    } while(inc >= 1);
    
    
    
    if(maxV>=10000) System.out.println("you should buy the absolute biggest reactor you can find");
    System.out.println("the best volume is: " +maxV+ " for the production of "+prodspec+". It creates a flow of "+maxFlow);

} 
                             
                             
                           
if(option=='S'){
   int tonum1=0;
   int tonum2=0;
   boolean flagger=true;
   char numerator='z';
   char denominator='z';
   while(flagger){
   denominator=safeGetChar("in terms of selectivity, input the species to be the denominator (it has to be a product)");
   tonum1= (int) denominator-'A';
   for(int z=0; z<reactstoic.length; z++){
     if(reactstoic[z][tonum1]>0) flagger=false;}
   }
   flagger=true;
   while(flagger){
   numerator=safeGetChar("in terms of selectivity, input the species to be the numerator (it has to be a product)");
   tonum2= (int) numerator-'A';
   for(int z=0; z<reactstoic.length; z++){
     if(reactstoic[z][tonum2]>0) flagger=false;}
   }

     double maxSel;
        double nextSel;
        double inc = 10000;
        double cap = 100000;
        double loBound = 0;
        double hiBound = cap;
        double maxV = 0;
        int initialMaxV = 0;
                
        do{
          maxSel = -1;
          for(v = loBound; v <= hiBound; v = v + inc)
          {
            rateInTermsFlow rate= new rateInTermsFlow();
            PFR[] pfr = new PFR[flow.length];
            
            for(int i=0; i<flow.length; i++) {
              pfr[i]=new PFR(i);
              pfr[i].setVals(flow, reactstoic, forwardk, keq, base, ct0, numrxn, rate, v, kd, rxnorder, phase, initflows);
            }
            SystemOfODEs client = new SystemOfODEs(pfr);
            double[] flow2=client.solveFlowRates(flow, v);

            
            nextSel = flow2[tonum1]/flow2[tonum2];
            
            if(nextSel > maxSel && Math.abs(nextSel-maxSel)/inc > 0.01)
            {
              maxSel = nextSel;
              maxV = v;
              // System.out.println("at: "+v+" "+"max flow: "+maxFlow);
            }
            else
              break;
          }
          if(maxV==0)
            loBound=10;
          else
            loBound = maxV - inc;
          if(maxV==cap)
            hiBound=cap;
          else
            hiBound = maxV + inc;
          inc = inc/2;
        } while(inc >= 1);
    
    
        if(maxV>=10000) System.out.println("you should buy the absolute biggest reactor you can find");
        else if((int)maxV == initialMaxV) System.out.println("The reactor volume should be minimized to maximize selectivity.");
        else System.out.println("the best volume is: " +maxV+ " for maximum selectivity of C/E, which is "+maxSel);
     }
    }
                                  
  }
}