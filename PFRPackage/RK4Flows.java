package PFRPackage;
import java.util.*;

public class RK4Flows{

	public double[] iterFlows(boolean isLiq,double[] flows,double[][] spec_coeffs,int rxNum,double[][] rxnK,double[] kInTermsOf,double[] interstage,int pfr_total,double volume){
		/*Uses adaptive RK4 to compute outlet molar flow rates*/
		flows[1] = interstage[0];
		double tot_flow = 0;
		for(int i = 0;i<flows.length;i++){
			tot_flow += flows[i];
		}
		Species spec_object = new Species(0.8,tot_flow,isLiq,flows,spec_coeffs);
		double[] fp1 = new double[spec_coeffs.length];

		double counter = 0;
		int pfr_num = 0;
		while(pfr_num<pfr_total){
			double[] r = spec_object.calcRi(spec_object.calculateConcs(spec_object.calcTotalFlow()),rxNum,rxnK,kInTermsOf);

			RK4_PFR rk4 = new RK4_PFR(spec_object,0.5);
			double[] k1 = rk4.calcK1(r);
			double[] k2 = rk4.calcK2calcK3(k1,rxnK,kInTermsOf,false);
			double[] k3 = rk4.calcK2calcK3(k2,rxnK,kInTermsOf,false);
			double[] k4 = rk4.calcK2calcK3(k3,rxnK,kInTermsOf,true);

			double[][] x = new double[flows.length][flows.length];
			for(int i = 0;i<x.length;i++){
				x[0][i] = k1[i];
				x[1][i] = k2[i];
				x[2][i] = k3[i];
				x[3][i] = k4[i];
			}//Hard-coded because number of rk4 k values are fixed.

			double[] y = rk4.nextYValues(x);

			fp1 = rk4.nextFlows(y);

			counter += 0.5;
			if(counter%volume == 0 && counter != 0){
				pfr_num +=1;
				if(pfr_num<pfr_total){
					fp1[1] += interstage[pfr_num];
				}
			}
			spec_object.setFlows(fp1);
		}

		return fp1;
	}

}