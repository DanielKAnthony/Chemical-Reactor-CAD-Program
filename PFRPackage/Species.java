package PFRPackage;
public class Species{
	private double c_T0;
	private double tot_flow;
	private boolean isLiq;
	private double[] flows;
	private double[][] coeffs;

	public Species(double c_T0, double tot_flow, boolean isLiq, double[] flows, double[][] coeffs) throws NullPointerException{
		try{
			this.c_T0 = c_T0;
			this.tot_flow = tot_flow;
			this.isLiq = isLiq;

			this.flows = new double[flows.length];
			for(int i = 0;i<flows.length;i++){
				this.flows[i] = flows[i];
			}

			this.coeffs = new double[coeffs.length][coeffs[0].length];
			for(int i = 0;i<coeffs.length;i++){
				for(int j = 0;j<coeffs[i].length;j++){
					this.coeffs[i][j] = coeffs[i][j];
				}
			}
		} catch(NullPointerException e){
			System.out.println("Error: Null value(s) present in data set");
		} catch(Exception e){
			System.out.println("An error occurred. Species object could not be constructed.");
		}
		
	}

	public Species(Species source){
		this.flows = new double[source.flows.length];
		for(int i = 0;i<source.flows.length;i++){
			this.flows[i] = source.flows[i];
		}

		this.coeffs = new double[source.coeffs.length][source.coeffs[0].length];
		for (int i = 0; i<source.coeffs.length;i++){
			for(int j = 0;j<source.coeffs[i].length;j++){
				this.coeffs[i][j] = source.coeffs[i][j];
			}
		}

		this.c_T0 = source.c_T0;
		this.tot_flow = source.tot_flow;
		this.isLiq = source.isLiq;
		
	}

	public Species clone(){
		return new Species(this);
	}

	public double[] getFlows(){
		double[] flows = new double[this.flows.length];
		for(int i = 0;i<flows.length;i++){
			flows[i] = this.flows[i];
		}
		return flows;
	}

	public double[][] getCoeffs(){
		double[][] coeffs = new double[this.coeffs.length][this.coeffs[0].length];
		for(int i = 0;i<this.coeffs.length;i++){
			for(int j =0;j<this.coeffs[i].length;j++){
				coeffs[i][j] = this.coeffs[i][j];
			}
		}
		return coeffs;
	}

	public boolean setFlows(double[] flows){
		if(flows.length != this.flows.length) return false;
		else if(flows == null) return false;
		for(int i = 0;i<flows.length;i++){
			this.flows[i] = flows[i];
		}
		return true;
	}

	public boolean setCoeffs(double[][] coeffs){
		if(coeffs.length != this.coeffs.length) return false;
		else if(coeffs == null) return false;
		for(int i = 0;i<coeffs.length;i++){
			if(coeffs[i] != this.coeffs[i]) return false;
			this.coeffs[i] = coeffs[i];
		}
		return true;
	}

	public boolean setCT0(double ct){
		this.c_T0 = ct;
		return true;
	}

	public double calcTotalFlow(){
		/*Calculates total flow. Updated if gas, if isLiq is true,
		returns fixed total flow.*/
		if(this.isLiq) return this.tot_flow;
		double flow_sum = 0;
		for (int i = 0;i<this.flows.length;i++){
			flow_sum += this.flows[i];
		}
		return flow_sum;
	}

	public double[] calculateConcs(double flow_sum){
		/*return array of concentrations for each species*/
		double[] speciesConcs = new double[this.flows.length];
		for (int i = 0;i<this.flows.length;i++){
			speciesConcs[i] = this.c_T0*(this.flows[i]/flow_sum);
		}
		return speciesConcs;
	}

		public double[][] getAllKValues(double[][] k, double[] coeffInTermsOf){
		/*Fills in k values for all other specs per reaction*/
		double[][] allKValues = new double[k.length][k[0].length];
		double[] activeK = new double[this.coeffs.length];

		for(int i = 0;i<k.length;i++){
			for(int j = 0;j<k[i].length;j++){
				if(k[i][j] != 0) activeK[i] = k[i][j];
			}
		}

		for(int i = 0;i<k.length;i++){
			for(int j = 0;j<k[i].length;j++){
				if(this.coeffs[i][j] != 0 && k[i][j] == 0){
					allKValues[i][j] = activeK[i]/coeffInTermsOf[i];
				}
				else if(this.coeffs[i][j] == 0){
					allKValues[i][j] = 0;
				}
				else if(k[i][j] != 0){
					allKValues[i][j] = k[i][j];
				}
			}
		}
		return allKValues;
	}

	public int countReactants(){
		/*Count number of reactants present.*/
		int reactant_count = 0;
		for(int i = 0;i<this.coeffs.length;i++){
			int temp_counter = 0;
			for(int j = 0;j<this.coeffs[i].length;j++){
				if((i == 0) && (this.coeffs[i][j]>0)) reactant_count += 1;
				else if(this.coeffs[i][j]>0) temp_counter += 1;
			}
			if (temp_counter>reactant_count) reactant_count = temp_counter;
		}
		return reactant_count;
	}

	public double[] IsolateReactantConcs(double[] concs){
		/*Ammend the returned concentration array of each species
		to solely represent reactant concentration. This reduces for-loop
		complexity when computing reaction rates.*/
		double[] reactant_concs = new double[concs.length];
		int reactant_count = countReactants();
		for(int i = 0;i<reactant_count;i++){
			reactant_concs[i] = concs[i];
		}
		for(int i = reactant_count;i<concs.length;i++){
			reactant_concs[i] = 0;
		}
		return reactant_concs;

	}

	public double[][] transposeMatrix(double[][] allKValues, double[][] k){
		/*Method used to get reaction constant matrices in sub-arrays that correspond
		to each species rather than by reaction*/
		double[][] temp_matrix = new double[k[0].length][k.length];
		for(int i = 0;i<k.length;i++){
			for(int j = 0;j<k[0].length;j++){
				temp_matrix[j][i] = allKValues[i][j];
			}
		}
		return temp_matrix;
	}

	public double[][] transposeMatrix(){
		/*Overloaded method used to get stoich coefficient sub-arrays
		categorized by species rather than by reaction.*/
		double[][] temp_matrix = new double[this.coeffs[0].length][this.coeffs.length];
		for(int i = 0;i<this.coeffs.length;i++){
			for(int j = 0;j<this.coeffs[i].length;j++){
				temp_matrix[j][i] = Math.abs(this.coeffs[i][j]);
			}
		}
		return temp_matrix;
	}

	public int specCount(int j_index){
		/*Count number of reactions each species is present in.
		Executed in calcRi method.*/
		int counter = 0;
		for(int i = 0;i<this.coeffs.length;i++){
			if(this.coeffs[i][j_index] != 0) counter += 1;
		}
		return counter;
	}

	public double[] combineExpressions(double[] arr, double[] reactant_concs,double[][] rxnCoeffs){
		/*multiplies reaction rate array elements by concentrations and their
		respective exponents*/
		double[] expressions = new double[arr.length];
		for(int i = 0;i<arr.length;i++){
			expressions[i] = 1;
			for(int j = 0;j<arr.length;j++){
				if(reactant_concs[j] != 0){
					expressions[i] *= Math.pow(reactant_concs[j],rxnCoeffs[i][j]);
				} 
			}
			arr[i] *= expressions[i];
		}
		return arr;
	}

	public double[] calcRi(double[] concs, int rxns, double[][] k, double[] coeffInTermsOf){
		/*returns array of reaction rates for each species.*/
		//Get all Ks for all reactions
		double[][] allKValues = getAllKValues(k,coeffInTermsOf);
		int reactant_count = countReactants();
		double[] reactant_concs = IsolateReactantConcs(concs);

		double[][] specKVals = transposeMatrix(allKValues,k);//arrays of k values spec by spec
		double[][] specCoeffs = transposeMatrix();//arrays of coeffs spec by spec

		
		//Repeat rx coefficients
		double[][] rxnCoeffs = new double[specCoeffs.length][specCoeffs[0].length];
		int counter = 0;

		for(int i = 0;i<reactant_count;i++){
			for(int j = 0;j<specCoeffs[i].length;j++){
				rxnCoeffs[i][j] = specCoeffs[i][j];
			}
		}

		for(int i = reactant_count; i<specCoeffs.length;i++){
			for(int j = 0;j<specCoeffs[i].length;j++){
				rxnCoeffs[i][j] = specCoeffs[0+counter][j];
			}
			counter += 1;

		}

		int[] speciesCounters = new int[concs.length];

		for(int i = 0;i<speciesCounters.length;i++){
			speciesCounters[i] = specCount(i);
		}

		double[] r_i = new double[specKVals.length];
		double[] temp = new double[r_i.length];
		int temp_index = 0;

		for(int i = 0;i<specKVals.length;i++){
			r_i[i] = 0;
			for(int j = 0;j<specKVals[i].length;j++){
				if(speciesCounters[i]<= 1){
					r_i[i] += specKVals[i][j];
				}
			
				else{
					temp[temp_index] = specKVals[i][j];
					temp_index+=1;
				}
			}
		}

		r_i = combineExpressions(r_i,reactant_concs,rxnCoeffs);
		temp = combineExpressions(temp,reactant_concs,rxnCoeffs);

		int index_calibrate = 0;
		for(int i = 0;i<reactant_count;i++){
			for (int j = 0;j<reactant_count;j++) {
				r_i[i] += temp[0+index_calibrate];
				index_calibrate += 1;	
			}
		}

		for(int i = 0;i<reactant_count;i++){
			r_i[i] = -1*r_i[i];
		}

		return r_i;
	}

}