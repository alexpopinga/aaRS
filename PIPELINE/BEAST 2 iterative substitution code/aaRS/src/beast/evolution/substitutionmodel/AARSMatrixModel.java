package beast.evolution.substitutionmodel;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;

import beast.core.Function;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;

public class AARSMatrixModel extends TriStateAARSSubstitutionModel implements Function {
	double distance;
	
	int [][][] m_R;
	
	@Override
	public void initAndValidate() {
		super.initAndValidate();
		m_R = new int[nrOfStates][][];
		for (int i = 0; i < nrOfStates; i++) {
			m_R[i] = new int[nrOfStates][nrOfStates];
			for (int j = 0; j < nrOfStates; j++) {
				Arrays.fill(m_R[i][j], -1);
			}
		}
	}
	
	
	@Override
	public int getDimension() {
		return 1;
	}

	@Override
	public double getArrayValue() {
		if (m_bRecalc) {
			update();
		}
		return distance;
	}

	@Override
	public double getArrayValue(int iDim) {
		if (m_bRecalc) {
			update();
		}
		return distance;
	}


	void update() {
		Tree tree = m_tree.get();
		// calculate epochs
		//Node[] nodes = tree.getNodesAsArray();
		
		// at this point, V (= internalnodes) contains nodes sorted by height
		V = new ArrayList<Node>();
		// start with any of the leaf nodes
		V.add(tree.getNode(0));

		calcStates(tree.getRoot(), V);

		Node thresholdNode = new Node();
		thresholdNode.setNr(1); // ensure we use m_fP[1] not used by other nodes
		double threshold = calcThreshold(tree.getRoot());
		thresholdNode.setHeight(threshold);
		V.add(thresholdNode);

		Collections.sort(V, comparator);
		
		// sanity check
		if (V.get(0) != tree.getNode(0) ||
				V.get(1) != thresholdNode) {
			System.err.println("Nodes in V are not as expected");
		}

		
		// initialise
		double[] fNewQ = m_fQ[0];		int [][] rNew = m_R[0];
		rNew[1][0] = m_nStateOfEpoch[tree.getRoot().getNr()];
		rNew[0][1] = rNew[1][0];

//		Arrays.fill(fNewQ, -1);
//		fNewQ[1 * nrOfStates + 0] = m_nStateOfEpoch[tree.getRoot().getNr()];
//		fNewQ[0 * nrOfStates + 1] = m_nStateOfEpoch[tree.getRoot().getNr()];

//		for n=2,19 % this is the dimension of the existing matrix
		for (int n2 = 2; n2 < V.size(); n2++) {
			int n = n2;
//			for i=n+1:I(n+1):-1
			int [][] rOld = rNew;
			rNew = m_R[n-1];
					
//			double [] fOldQ = fNewQ;
//			fNewQ = m_fQ[n];
			Node node = V.get(n);
			int iStateLeft = m_nStateOfEpoch[node.getLeft().getNr()];
			int iStateRight = m_nStateOfEpoch[node.getRight().getNr()];
			int i2 = Math.min(iStateLeft, iStateRight);
			int i1 = Math.max(iStateLeft, iStateRight);
			
			
			// copy
			for (int i = 0; i < nrOfStates; i++) {
				for (int j = 0; j < nrOfStates; j++) {
					rNew[i][j] = rOld[i][j];
				}
			}
			// i1 row
			for (int j = 0; j < nrOfStates; j++) {
				rNew[i1][j] = rOld[i2][j];
			}
			// i1 col
			for (int i = 0; i < nrOfStates; i++) {
				rNew[i][i1] = rOld[i][i2];
			}
			// new entry
			rNew[i1][i1] = -1;
			rNew[i2][i1] = i1;
			rNew[i1][i2] = i1;

			
//			// upper left corner
//			for (int i = 0; i < i1; i++) {
//				for (int j = 0; j < i1; j++) {
//					rNew[i][j] = rOld[i][j];
//					
//				}
//			}
//			// upper right corner
//			for (int i = i1; i < n; i++) {
//				for (int j = 0; j < i1; j++) {
//					rNew[i][j] = rOld[i - 1][j];
//					
//				}
//			}
//			// lower left corner
//			for (int i = 0; i < i1; i++) {
//				for (int j = i1; j < n; j++) {
//					rNew[i][j] = rOld[i][j-1];
//					
//				}
//			}
//			// lower right corner
//			for (int i = i1; i < n; i++) {
//				for (int j = i1; j < n; j++) {
//					rNew[i][j] = rOld[i - 1][j-1];
//					
//				}
//			}
//			
//			// iMin row
//			for (int j = 0; j < i1; j++) {rNew[i1][j] = rOld[i2][j];}
//			rNew[i1][i1] = 0;
//			for (int j = i1; j < n; j++) {rNew[i1][j] = rOld[i2][j-1];}
//			
//			// iMin column
//			for (int i = 0; i < i1; i++) {rNew[i][i1] = rOld[i][i2];}
//			rNew[i1][i1] = 0;
//			for (int i = i1; i < n; i++) {rNew[i][i1] = rOld[i-1][i2];}
//			
//			// new entry
//			rNew[i2][i1] = i1;
//			rNew[i1][i2] = i1;
			

//			for (int i = 0; i < iMax; i++) {
//				System.arraycopy(fOldQ, nrOfStates * i, fNewQ, nrOfStates * i, iMax);
//				fNewQ[nrOfStates * i + iMax] = fOldQ[nrOfStates * i + iMin];
//				System.arraycopy(fOldQ, nrOfStates * i + iMax, fNewQ, nrOfStates * i + iMax + 1, nrOfStates - iMax - 1);
//			}
//			System.arraycopy(fOldQ, nrOfStates * iMin, fNewQ, nrOfStates * iMax, iMax);
//			fNewQ[nrOfStates * iMax + iMax] = -1;
//			System.arraycopy(fOldQ, nrOfStates * iMin, fNewQ, nrOfStates * iMax + 1, nrOfStates - iMax - 1);
//			for (int i = iMax + 1; i < nrOfStates; i++) {
//				System.arraycopy(fOldQ, nrOfStates * (i - 1), fNewQ, nrOfStates * i, iMax);
//				fNewQ[nrOfStates * (i + 1) + iMax] = fOldQ[nrOfStates * i + iMin];
//				System.arraycopy(fOldQ, nrOfStates * (i - 1) + iMax, fNewQ, nrOfStates * i + iMax + 1, nrOfStates - iMax - 1);
//			}
//			fNewQ[nrOfStates * iMin + iMax] = iMax;
//			fNewQ[nrOfStates * iMax + iMin] = iMax;
			
//			for (int i = n + 1; i >= ; i--) {
//				A(i)=A(i-1)
//				for j=1:n
//					% this copies the RHS of the matrix across
//					% i is the index for a column
//					% j is the index for a row
//					R(i,j)=R(i-1,j)
//				end
//			end
//			for j=n+1:I(n):-1 %now copy the bottom right quadrant down
//				for i=n+1:I(n+1):-1
//				% i is the index for a row; j for column
//				R(i,j)=R(i,j-1)
//				end
//			end
//			R(I(k),I(k+1))=alpha(k+1)
//			% this occupies the new “blank” space
//		end	
		}				  

		
		fit = calcFit();
		m_bRecalc = false;
	}

	public double calcFit() {
		SubstitutionModel model = m_default.get();
		double [][] matrix = ((GeneralSubstitutionModel)model).getRateMatrix();
		double [] weight = new double[nrOfStates];
		double [] sumDiff = new double[nrOfStates];
		int [][] indexMatrix = m_R[nrOfStates - 1]; 
		for (int i = 0; i < nrOfStates * nrOfStates; i++) {
			int index = indexMatrix[i/nrOfStates][i%nrOfStates];
			if (index != -1) {
				weight[index] += 1;
				sumDiff[index] += matrix[i/nrOfStates][i%nrOfStates];
			}
		}
		double [] alpha = new double[nrOfStates];
		for (int i = 0; i < nrOfStates; i++) {
			alpha[i] = sumDiff[i] / weight[i]; 
		}

		double misfit = 0;
		for (int i = 0; i < nrOfStates * nrOfStates; i++) {
			int index = indexMatrix[i/nrOfStates][i%nrOfStates];;
			if (index != -1) {
				double rate = matrix[i/nrOfStates][i%nrOfStates];
				misfit += (alpha[index] - rate) * (alpha[index] -rate); 
			}
		}
		return misfit;
	}

	double fit;
	public double getFit() {
		if (m_bRecalc) {
			update();
		}
		return fit;
	}
}
