package beast.evolution.substitutionmodel;

import test.beast.beast2vs1.TreePriorTest;
import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.evolution.tree.Node;
import beast.evolution.tree.TraitSet;
import beast.evolution.tree.Tree;
import beast.evolution.tree.TreeDistribution;
import beast.evolution.tree.coalescent.TreeIntervals;

@Description("Prior ensuring all organism trees are below the AARS tree")
public class AARSTreePrior extends TreeDistribution {
	public Input<TraitSet> m_traits = new Input<TraitSet>("trait", "traiset describing which taxon goes with which AARS", Validate.REQUIRED); 

	TreeIntervals m_intervals;
	int [] m_nStateOfEpoch;
	
	@Override
	public void initAndValidate() {
		m_intervals = treeIntervalsInput.get();
		if (m_intervals == null) {
			throw new RuntimeException("expected treeintervals to be specified");
		}
		m_nStateOfEpoch = new int[m_intervals.treeInput.get().getNodeCount()];
	}
	
	final private double THRESHOLD = 1;
	final private double ZERO = 1e-6;
	


	@Override
	public double calculateLogP() {
		// ensure AARS nodes are above organism nodes
		Tree tree = m_intervals.treeInput.get();
		calcStates(tree.getRoot());
		Node [] nodes = tree.getNodesAsArray();
		double fMaxOrganismHeight = 0;
		double fMinAARSHeight = Double.MAX_VALUE;
		for (int i = 0; i < nodes.length; i++) {
			if (m_nStateOfEpoch[i] < 0) {
				fMinAARSHeight = Math.min(fMinAARSHeight, nodes[i].getHeight());
			} else {
				fMaxOrganismHeight = Math.max(fMaxOrganismHeight, nodes[i].getHeight());
			}
		}
		if (fMaxOrganismHeight > fMinAARSHeight) {
			return Double.NEGATIVE_INFINITY;
		}
		
		
		logP = 0;
		if (true) {
			return 0;
		}
		
		// now, ensure AARS nodes don't clump 
		double [] inters = m_intervals.getIntervals(null);
		//for (double diff :  m_intervals.getIntervals(null)) {
//			if (diff > 0 && diff < 1) {
//				logP += /* 5* */ Math.log(diff);
//			}
		for (int i = 0; i < 20; i++) {
			double diff = inters[inters.length-1-i];
			if (diff <= ZERO) {
				logP += 100 * Math.log(ZERO/THRESHOLD);
			} else if (diff < THRESHOLD) {
				logP += 100 * Math.log(diff/THRESHOLD);
			}
		}
		return logP;
	}
	
	int calcStates(Node node) {
		//m_fDepths[node.getNr()] = node.getHeight();
		if (node.isLeaf()) {
			m_nStateOfEpoch[node.getNr()] = (int) m_traits.get().getValue(node.getNr());
		} else {
			int iStateLeft = calcStates(node.getLeft());
			int iStateRight = calcStates(node.getRight());
			m_nStateOfEpoch[node.getNr()] = Math.min(iStateLeft, iStateRight);
			if (iStateLeft != iStateRight) {
				m_nStateOfEpoch[node.getNr()] = -1;
			}
		}
		return m_nStateOfEpoch[node.getNr()];
	}
	
}
