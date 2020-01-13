package PFRPackage;
import java.util.*;

public class FlowSolver{

	public int factorial(int n){
	    if (n == 0) return 1;
	                
	    return n*factorial(n-1); 
	}//factorial of number of reactors using recursion


    public double[][] getAllPermutations(double[] flow_pattern){
    	/*permute every new flow configuration that is generated in the
    	flows method*/
	    List<Integer> current = new ArrayList<Integer>();
	    List<Integer[]> permutations = new ArrayList<Integer[]>();

	    int length = flow_pattern.length;
	    current.add(-1);

	    while (!current.isEmpty()) {
	        int position = Integer.MAX_VALUE;
	        position = nextAvailable(current, current.get(current.size()-1) + 1);
	        current.remove(current.size()-1);
	        while (position >= length && !current.isEmpty()) {
	            position = nextAvailable(current, current.get(current.size()-1) + 1);
	            current.remove(current.size()-1);
	        }
	        if (position < length) {
	            current.add(position);
	        } else {
	            break;
	        }

	        while (current.size() < length) {
	            int unused = nextAvailable(current, 0);

	            current.add(unused);
	        }
	        permutations.add(current.toArray(new Integer[0]));
	    }

	    int numPermutations = permutations.size();
	    double[][] result = new double[numPermutations][length];
	    for (int i = 0; i < numPermutations; ++i) {
	        Integer[] indexes = permutations.get(i);
	        double[] row = new double[length];
	        for (int d = 0;d<length; ++d) {
	            row[d] = flow_pattern[indexes[d]];
	        }
	        result[i] = row;
	    }

	    return result;
	}//Returns array containing arrays of of all permutations per flow config.

	public int nextAvailable(List<Integer> used, Integer current){
		/*gets next array element that is not present in current permuted sequence.
		This is a helper method for getAllPermutations*/
	    int unused = current != null ? current : 0;
	    while (used.contains(unused)) {
	        ++unused;
	    }
	    return unused;
	}//gets next array element that is not present in current permuted sequence

	 public double computeError(int max_flow){
	  double error = max_flow/2;
	  boolean flag = true;
	  while(error>0.005*max_flow){
	    error -= 0.01;
	  }
	  return error;

	 }//computes error tolerance for next flow iteration

	 public double[][] flows(int max_flow,int rxrs){
	  /*returns 2d array of all permuted flow pattern for interstage feeds*/
	  double[][] configs = new double[((int)(max_flow/computeError(max_flow))+1)*factorial(rxrs)][rxrs];
	  double[] b = new double[rxrs];
	  b[0] = max_flow;
	  for(int i = 1;i<b.length;i++){
	   b[i] = 0;
	  }
	  double[] init = Arrays.copyOf(b,b.length);
	  configs[0] = init;

	  int i = 1;
	  boolean flag = true;
	  int c_index = 1;
	  while(flag){
	   b[0] = Math.round((b[0]-computeError(max_flow))*10000.0)/10000.0;
	   b[i] = Math.round((b[i]+computeError(max_flow))*10000.0)/10000.0;
	   if((i+1) == b.length){
	    i = 1;
	   }
	   else{
	    i += 1;
	   }
	   double flow_sum = 0;
	   
	   for(int j = 0;j<b.length;j++){
	    flow_sum += b[j];
	   }
	   if((int)Math.round(flow_sum) == (int)max_flow){
	    double[] con = Arrays.copyOf(b,b.length);
	    double[][] temp_arr = getAllPermutations(con);
	    for(int k = 0;k<temp_arr.length;k++){
	     configs[c_index] = temp_arr[k];
	     c_index += 1;
	    }

	   }

	   if(b[0] == 0){
	    flag = false;
	   }

	  }
	  return configs;

	 }
}