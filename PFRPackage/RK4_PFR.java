package PFRPackage;
import java.util.*;
import java.io.*;

public class RK4_PFR extends Species{

	private double step_size;

	public RK4_PFR(Species specParams, double step_size){
		super(specParams);
		this.step_size = step_size;
	}

	public RK4_PFR(RK4_PFR source){
		super(source);
		this.step_size = source.step_size;
	}

	public RK4_PFR clone(){
		return new RK4_PFR(this);
	}

	public double getStepSize(){
		return this.step_size;
	}

	public boolean setStepSize(double step_size){
		this.step_size = step_size;
		return true;
	}

	public double[] combineExps(double[] arr, double[] rxnFlows,double[] rk4_ki,int flow_index,int exp_count){
		/*Combines rk4 expressions to yield each k2,k3 term*/
		double[] expressions = new double[arr.length];
		for(int i = 0;i<arr.length-1;i++){
			expressions[i] = 1;
			if(exp_count == 1){
				expressions[i] *= Math.pow((rxnFlows[flow_index]+((this.step_size/2)*rk4_ki[flow_index]))/super.calcTotalFlow(),2);
				expressions[i] *= (rxnFlows[flow_index+1]+((this.step_size/2)*rk4_ki[flow_index+1]))/super.calcTotalFlow();
			}
			else{
				expressions[i] *= (rxnFlows[flow_index]+((this.step_size/2)*rk4_ki[flow_index]))/super.calcTotalFlow();
				expressions[i] *= Math.pow((rxnFlows[flow_index+1]+((this.step_size/2)*rk4_ki[flow_index+1]))/super.calcTotalFlow(),2);
			}
			arr[i] *= expressions[i];
		}

		arr[arr.length-1] *= (rxnFlows[flow_index]+((this.step_size/2)*rk4_ki[flow_index]))/super.calcTotalFlow();
		arr[arr.length-1] *= Math.pow((rxnFlows[flow_index+1]+((this.step_size/2)*rk4_ki[flow_index+1]))/super.calcTotalFlow(),2);

		return arr;
	}

	public double[] combineExpsK4(double[] arr, double[] rxnFlows,double[] rk4_ki,int flow_index,int exp_count){
		/*Same for k4 term*/
		double[] expressions = new double[arr.length];
		for(int i = 0;i<arr.length-1;i++){
			expressions[i] = 1;
			if(exp_count == 1){
				expressions[i] *= Math.pow((rxnFlows[flow_index]+rk4_ki[flow_index])/super.calcTotalFlow(),2);
				expressions[i] *= (rxnFlows[flow_index+1]+rk4_ki[flow_index+1])/super.calcTotalFlow();
			}
			else{
				expressions[i] *= (rxnFlows[flow_index]+rk4_ki[flow_index])/super.calcTotalFlow();
				expressions[i] *= Math.pow((rxnFlows[flow_index+1]+rk4_ki[flow_index+1])/super.calcTotalFlow(),2);
			}
			arr[i] *= expressions[i];
		}

		arr[arr.length-1] *= (rxnFlows[flow_index]+rk4_ki[flow_index])/super.calcTotalFlow();
		arr[arr.length-1] *= Math.pow((rxnFlows[flow_index+1]+rk4_ki[flow_index+1])/super.calcTotalFlow(),2);

		return arr;
	}

	public double[] calcK1(double[] r_i){
		/*Get k1 terms*/
		return r_i;
	}

	public double[] calcK2calcK3(double[] rk4_ki,double[][] k, double[] coeffInTermsOf, boolean isK4){
		/*Get k2, k3, k4 terms*/
		double[][] allKValues = getAllKValues(k,coeffInTermsOf);
		double[][] specKVals = transposeMatrix(allKValues,k);//arrays of k values spec by spec


		double[] ki_minus1 = new double[rk4_ki.length];
		double[] reactant_flows = new double[rk4_ki.length];

		for(int i = 0;i<super.countReactants();i++){
			ki_minus1[i] = rk4_ki[i];	
			reactant_flows[i] = super.getFlows()[i];
		}//isolate previous k values and flows for reactants.

		int inc = 0;

		for(int i = ki_minus1.length;i<rk4_ki.length;i++){
			ki_minus1[i] = ki_minus1[0+inc];
			reactant_flows[i] = reactant_flows[0+inc];
			inc += 1;
		}

		int[] speciesCounters = new int[rk4_ki.length];

		for(int i = 0;i<speciesCounters.length;i++){
			speciesCounters[i] = specCount(i);
		}

		double[] k23i = new double[rk4_ki.length];
		double[] temp = new double[k23i.length];
		int temp_index = 0;

		for(int i = 0;i<rk4_ki.length;i++){
			for(int j = 0;j<specKVals[i].length;j++){
				if(speciesCounters[i] <= 1){
					k23i[i] += specKVals[i][j];
				}
				else{
					temp[i] = specKVals[i][j];
					temp_index += 1;		
				}
			}

		}

		for(int i = 0;i<super.countReactants();i++){
				k23i[i] = specKVals[i][0];
		}

		if(!isK4){
			k23i = combineExps(k23i,reactant_flows,ki_minus1,0,1);
			temp = combineExps(temp,reactant_flows,ki_minus1,0,2);
		}

		else{

			k23i = combineExpsK4(k23i,reactant_flows,ki_minus1,0,1);
			temp = combineExpsK4(temp,reactant_flows,ki_minus1,0,2);

		}

		double[] multipleExpressions = new double[super.countReactants()];
		for(int i = 0;i<multipleExpressions.length;i++){
			multipleExpressions[i] = k23i[i] + temp[i];
			k23i[i] = multipleExpressions[i];
		}

		for(int i = 0;i<super.countReactants();i++){
			k23i[i] *= -1;
		}

		return k23i;

	}

	public double[] nextYValues(double[][] rk4_k){
		double[] y_iPlus1 = new double[super.getFlows().length];

		for(int i = 0; i<y_iPlus1.length;i++){
			y_iPlus1[i] = 0;
			for(int j = 0;j<rk4_k[0].length;j++){
				if(j == 1 || j == 2) y_iPlus1[i] += 2*rk4_k[j][i];
				else y_iPlus1[i] += rk4_k[j][i];
			}
			y_iPlus1[i] /= 6;
		}

		return y_iPlus1;
	}

	public double[] nextFlows(double[] y_iPlus1){
		double[] fPlus1 = new double[y_iPlus1.length];
		for(int i = 0;i<y_iPlus1.length;i++){
			fPlus1[i] = super.getFlows()[i]+(y_iPlus1[i]*this.step_size);
		}

		return fPlus1;
	}


}