import java.io.*;
import java.util.*;
import RK4Package.*;
import PFRPackage.*;

public class FinalDriver{
 public static void main(String[] args) throws FileNotFoundException,InterruptedException {
  System.out.println("_______Group 14 Reactor Simulator_______\n");
  Scanner scan = new Scanner(System.in);
  boolean on_repeat = true;
  while(on_repeat){

   boolean choice_check = true;
   int user_choice = -1;

   while(choice_check){
    try{
     System.out.println("Enter 1 to run a pre-programmed solution, or enter 0 to manually enter a data set");
     user_choice = scan.nextInt();
     if(user_choice == 1 || user_choice == 0) choice_check = false;
     else System.out.println("Please enter either 1 or 0");
    }catch(Exception e){
     System.out.println("Invalid");
     scan.nextLine();
    }
   }

   switch(user_choice){
    
    case 1:
     boolean file_check = true;
     int file_choice = -1;
     while(file_check){
      try{
       System.out.println("Enter 0 to run the base case\nEnter 1 to run Problem 1\nEnter 2 to run Problem 2");
       file_choice = scan.nextInt();
       if(file_choice<=2 && file_choice>=0) file_check = false;
       else System.out.println("Input is out of range");
      }catch(Exception e){
       System.out.println("Invalid");
       scan.nextLine();
      }
     }

     if(file_choice == 0) autoFill(0);
     else if(file_choice == 1) autoFill(1);
     else autoFill(2);
    case 0:
     if(user_choice == 0) manualEntry();
     // manualEntry(); 
   }

   int rep = -1;
   boolean rep_check = true;
   while(rep_check){
    try{
     System.out.println("\nEnter 1 to perform another simulation\nEnter 0 to exit the program");
     rep = scan.nextInt();
     if(rep == 0){
      on_repeat = false;
      break;
     }
     else if(rep == 1) break;
     else System.out.println("Enter 1 or 0");
    }catch(Exception e){
     System.out.println("Invalid");
     scan.nextLine();
    }
   }
  }
 }

 public static void autoFill(int n) throws FileNotFoundException, InterruptedException {
  // FileInputStream ins = null;
  System.out.println("Retrieving data...");
  Thread.sleep(250);

  if(n == 0){
   FileInputStream pb0 = null;
   try{
    pb0 = new FileInputStream("Problem0_Group14.txt");
   }catch(FileNotFoundException e){
    System.out.println("The required file could not be loaded.");
   }
   Scanner inputs = new Scanner(pb0);
   double[] flows = new double[3];
   for(int i = 0;i<flows.length;i++){
    flows[i] = inputs.nextDouble();
   }
   double[][] stoich = new double[2][3];
   for(int i = 0;i<stoich.length;i++){
    for(int j = 0;j<stoich[i].length;j++){
     stoich[i][j] = inputs.nextDouble();
    }
   }
   int numofrxn = inputs.nextInt();
   double ct0 = inputs.nextDouble();
   double[] k = new double[2];
   for(int i = 0;i<k.length;i++){
    k[i] = inputs.nextDouble();
   }
   int[] base = new int[2];
   for(int i = 0;i<base.length;i++){
    base[i] = inputs.nextInt();
   }
   double v = inputs.nextDouble();
   double[] keq = new double[2];
   for(int i = 0;i<keq.length;i++){
    keq[i] = inputs.nextDouble();
   }

   double[] kd = new double[]{0,0,0};
   rateInTermsFlow rate= new rateInTermsFlow();
            PFR[] pfr = new PFR[flows.length];
            double[][] rxnorder = stoich;

            for(int i=0; i<flows.length; i++) {
               pfr[i]=new PFR(i);
               pfr[i].setVals(flows, stoich, k, keq, base, ct0, 1, rate, v, kd, rxnorder, true,flows);
            }
            SystemOfODEs client = new SystemOfODEs(pfr);


             double[] flow2=client.solveFlowRates(flows, v);
             char c = 'A';
             for(int i=0; i<flow2.length; i++){
               System.out.println( "flow is :" +flow2[i]+"  for species "+c);
               c++;
             }

  }
  else if(n == 1){
   FileInputStream pb1 = null;
   try{
    pb1 = new FileInputStream("Problem1_Group14.txt");
   }catch(FileNotFoundException e){
    System.out.println("The required file could not be loaded.");
   }
   Scanner inputs = new Scanner(pb1);

   int specs = inputs.nextInt();
   int rxns = inputs.nextInt();
   double[][] coeffs = new double[rxns][specs];
   for(int i = 0;i<coeffs.length;i++){
    for(int j = 0;j<coeffs[i].length;j++){
     coeffs[i][j] = inputs.nextDouble();
    }
   }
   int kinetics = inputs.nextInt();//elementary
   double[] keq = new double[rxns];
   for(int i = 0;i<keq.length;i++){
    keq[i] = inputs.nextDouble();
   }
   double[] rates = new double[keq.length];
   for(int i = 0;i<rates.length;i++){
    rates[i] = inputs.nextDouble();
   }
   double[] diffusivities = new double[specs];
   for(int i = 0;i<specs;i++){
    diffusivities[i] = inputs.nextDouble();
   }
   double[] flows = new double[specs];
   for(int i = 0;i<specs;i++){
    flows[i] = inputs.nextDouble();
   }
   double v_0 = inputs.nextDouble();

   double ct0 = 30/v_0;
   int[] base = new int[]{0,2};
   double[][] rxnorder = new double[coeffs.length][coeffs[0].length];
   for (int i=0; i<coeffs.length; i++){
      for(int j=0; j<coeffs[i].length; j++)
       rxnorder[i][j]=coeffs[i][j];
     }

   //Prod
    double maxFlow;
      double nextFlow;
      double inc = 10000;
      double cap = 100000;
      double loBound = 0;
      double hiBound = cap;
      double maxV = 0;
      int tonum = 2;
      double[] temp = new double[]{0,0};
      
      do{
        maxFlow = -1;
        for(double v = loBound; v <= hiBound; v = v + inc)
        {
          rateInTermsFlow rate= new rateInTermsFlow();
          PFR[] pfr = new PFR[flows.length];
          
          for(int i=0; i<flows.length; i++) {
            pfr[i]=new PFR(i);
            pfr[i].setVals(flows, coeffs, rates, keq, base, ct0, rxns, rate, v_0, diffusivities, rxnorder, false,temp);
          }
          SystemOfODEs client = new SystemOfODEs(pfr);
          double[] flow2=client.solveFlowRates(flows, v);
          
          nextFlow = flow2[tonum];
          
          if(nextFlow > maxFlow && Math.abs(nextFlow-maxFlow)/inc > 0.001)
          {
            maxFlow = nextFlow;
            maxV = v;
            // System.out.println("at: "+v+" "+"max flow: "+maxFlow);
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
      System.out.println("the best volume is: " +maxV+ " for the production of C. It creates a flow of "+maxFlow);

   //Sel
      double maxSel;
      double nextSel;
      inc = 10000;
      cap = 100000;
      loBound = 0;
      hiBound = cap;
      maxV = 0;
      int initialMaxV = 0;
      // int tonum = 2;
      int tonum1 = 2;
      int tonum2 = 4;
      
      do{
        maxSel = -1;
        for(double v = loBound; v <= hiBound; v = v + inc)
        {
          rateInTermsFlow rate= new rateInTermsFlow();
          PFR[] pfr = new PFR[flows.length];
          
          for(int i=0; i<flows.length; i++) {
            pfr[i]=new PFR(i);
            pfr[i].setVals(flows, coeffs, rates, keq, base, ct0, rxns, rate, v_0, diffusivities, rxnorder, false,temp);
          }
          SystemOfODEs client = new SystemOfODEs(pfr);
          double[] flow2=client.solveFlowRates(flows, v);

          
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


  else if(n == 2){
   
   FileInputStream pb2 = null;
   try{
    pb2 = new FileInputStream("Problem2_Group14.txt");
   }catch(FileNotFoundException e){
    System.out.println("The required file could not be loaded.");
   }
   Scanner inputs = new Scanner(pb2);
   int rxnum = inputs.nextInt();
   int specs = inputs.nextInt();
   double[][] spec_coeffs = new double[rxnum][specs];
   for(int i = 0;i<spec_coeffs.length;i++){
    for(int j = 0;j<spec_coeffs[i].length;j++){
     spec_coeffs[i][j] = inputs.nextDouble();
    }
   }
   double[] flows = new double[specs];
   for(int i = 0;i<flows.length;i++){
    flows[i] = inputs.nextDouble();
   }
   double[] userK = new double[rxnum];
   for(int i = 0;i<userK.length;i++){
    userK[i] = inputs.nextDouble();
   }
   int[] k_index = new int[rxnum];
   for(int i = 0;i<k_index.length;i++){
    k_index[i] = inputs.nextInt();
   }
   double[] kInTermsOf = new double[k_index.length];
   for(int i = 0;i<kInTermsOf.length;i++){
    kInTermsOf[i] = spec_coeffs[i][k_index[i]];
   }
   double[][] rxnK = new double[rxnum][specs];
   for(int i = 0;i<rxnum;i++){
    for(int j = 0;j<specs;j++){
     if(k_index[i] == j){
      rxnK[i][j] = userK[i];
     }
     else{
      rxnK[i][j] = 0;
     }
    }
   }

   inputs.close();

   RK4Flows rkflow = new RK4Flows();
   FlowSolver flow_solve = new FlowSolver();

   double[] interstage = new double[] {flows[1],1,1,1,1};
   double[] test = rkflow.iterFlows(false,flows,spec_coeffs,rxnum,rxnK,kInTermsOf,interstage,5,10);
   double[][] f = flow_solve.flows(4,5);

   List<Double> sels = new ArrayList<Double>();
   System.out.print("Optimizing selectivity using interstage feeds of species B\n\n45 seconds remaining\r");
   for(int i = 0;i<f.length;i++){
    sels.add(rkflow.iterFlows(false,flows,spec_coeffs,rxnum,rxnK,kInTermsOf,f[i],5,10)[2]/rkflow.iterFlows(false,flows,spec_coeffs,rxnum,rxnK,kInTermsOf,f[i],5,10)[3]);
   }

   int index_count = 0;
   double max = 0;
   for(int i = 0;i<sels.size();i++){
    if(sels.get(i)>max){
     max = sels.get(i);
     index_count = i;
    }
   }
   // System.out.print();
   System.out.print("Optimization complete\t\nPress ENTER to see optimized results");

   waitForUser();

   System.out.println("\nMaximum selectivity is " +max+ "\nusing Fb values of:");
   for(int i = 0;i<f[0].length;i++){
    System.out.print(f[index_count][i] + " mol/min into PFR "+(i+1)+"\n");
   }

   char react_char = 'A';
   System.out.println("\nThis optimized configuration produces the following outlet flows from the final PFR:\n");
   for(int i = 0;i<test.length;i++){
    System.out.println("F_"+react_char+" = "+rkflow.iterFlows(false,flows,spec_coeffs,rxnum,rxnK,kInTermsOf,f[index_count],5,10)[i]+" mol/min");
    react_char++;
   }
  }



 }

 public static void manualEntry(){
  Manual man = new Manual();
  man.main();
 }

 public static void waitForUser(){
  Scanner prompt = new Scanner(System.in);
  prompt.nextLine();
 }
}