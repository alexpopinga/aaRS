/*
* File SubtreeSlide.java
*
* Copyright (C) 2010 Remco Bouckaert remco@cs.auckland.ac.nz
*
* This file is part of BEAST2.
* See the NOTICE file distributed with this work for additional
* information regarding copyright ownership and licensing.
*
* BEAST is free software; you can redistribute it and/or modify
* it under the terms of the GNU Lesser General Public License as
* published by the Free Software Foundation; either version 2
* of the License, or (at your option) any later version.
*
*  BEAST is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU Lesser General Public License for more details.
*
* You should have received a copy of the GNU Lesser General Public
* License along with BEAST; if not, write to the
* Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
* Boston, MA  02110-1301  USA
*/
/*
 * SubtreeSlideOperator.java
 *
 * Copyright (C) 2002-2006 Alexei Drummond and Andrew Rambaut
 *
 * This file is part of BEAST.
 * See the NOTICE file distributed with this work for additional
 * information regarding copyright ownership and licensing.
 *
 * BEAST is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 *  BEAST is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 */

package beast.evolution.operators;


import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.evolution.substitutionmodel.AARSSubstitutionModel;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.util.Randomizer;

import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.List;

/**
 * Implements the subtree slide move.
 */
@Description("SubtreeSlide for AARS epoch nodes only. " +
		"Note that no attempt is made to restrict epoch nodes to slide among the other. " +
		"However, with the AARSConstraint, this should never happen. " +
		"Therefore, the Hasting ratio remain as for SubtreeSlide.")
public class AARSSubtreeSlide extends TreeOperator {

    public Input<Double> m_size = new Input<Double>("size", "size of the slide, default 1.0", 1.0);
    public Input<Boolean> m_gaussian = new Input<Boolean>("gaussian", "Gaussian (=true=default) or uniform delta", true);
    public Input<Boolean> m_bOptimise = new Input<Boolean>("optimise", "flag to indicate that the scale factor is automatically changed in order to acheive a good acceptance rate (default true)", true);

    public Input<AARSSubstitutionModel> substModelInput = new Input<AARSSubstitutionModel>("substModel","AARS substitution model providing list or Epoch nodes to operator on", Validate.REQUIRED);

    AARSSubstitutionModel substModel;
    
    // shadows m_size
    double m_fSize;

    @Override
    public void initAndValidate() {
    	substModel = substModelInput.get();
        m_fSize = m_size.get();
    }

    /**
     * Do a probabilistic subtree slide move.
     *
     * @return log of Hastings Ratio, or Double.NEGATIVE_INFINITY if proposal should not be accepted *
     */
    @Override
    public double proposal() {
        final Tree tree = treeInput.get(this);

        if (substModel.getV().size() < 19) {
        	return Double.NEGATIVE_INFINITY;
        }

        double logq;

        Node i;

        // 1. choose a random node avoiding root
        do {
            i = tree.getNode(substModel.getV().get(1+Randomizer.nextInt(substModel.getV().size()-1)).getNr());
            //i = tree.getNode(Randomizer.nextInt(tree.getNodeCount()));
        } while (i.isRoot());

        final Node iP = i.getParent();
        final Node CiP = getOtherChild(iP, i);
        final Node PiP = iP.getParent();

        // 2. choose a delta to move
        final double delta = getDelta();
        final double oldHeight = iP.getHeight();
        final double newHeight = oldHeight + delta;

        // 3. if the move is up
        if (delta > 0) {

            // 3.1 if the topology will change
            if (PiP != null && PiP.getHeight() < newHeight) {
                // find new parent
                Node newParent = PiP;
                Node newChild = iP;
                while (newParent.getHeight() < newHeight) {
                    newChild = newParent;
                    newParent = newParent.getParent();
                    if (newParent == null) break;
                }


                // 3.1.1 if creating a new root
                if (newChild.isRoot()) {
                    replace(iP, CiP, newChild);
                    replace(PiP, iP, CiP);

                    iP.setParent(null);
                    tree.setRoot(iP);

                }
                // 3.1.2 no new root
                else {
                    replace(iP, CiP, newChild);
                    replace(PiP, iP, CiP);
                    replace(newParent, newChild, iP);
                }

                iP.setHeight(newHeight);

                // 3.1.3 count the hypothetical sources of this destination.
                final int possibleSources = intersectingEdges(newChild, oldHeight, null);
                //System.out.println("possible sources = " + possibleSources);

                logq = -Math.log(possibleSources);

            } else {
                // just change the node height
                iP.setHeight(newHeight);
                logq = 0.0;
            }
        }
        // 4 if we are sliding the subtree down.
        else {

            // 4.0 is it a valid move?
            if (i.getHeight() > newHeight) {
                return Double.NEGATIVE_INFINITY;
            }

            // 4.1 will the move change the topology
            if (CiP.getHeight() > newHeight) {

                final List<Node> newChildren = new ArrayList<Node>();
                final int possibleDestinations = intersectingEdges(CiP, newHeight, newChildren);

                // if no valid destinations then return a failure
                if (newChildren.size() == 0) {
                    return Double.NEGATIVE_INFINITY;
                }

                // pick a random parent/child destination edge uniformly from options
                final int childIndex = Randomizer.nextInt(newChildren.size());
                final Node newChild = newChildren.get(childIndex);
                final Node newParent = newChild.getParent();


                // 4.1.1 if iP was root
                if (iP.isRoot()) {
                    // new root is CiP
                    replace(iP, CiP, newChild);
                    replace(newParent, newChild, iP);

                    CiP.setParent(null);
                    tree.setRoot(CiP);

                } else {
                    replace(iP, CiP, newChild);
                    replace(PiP, iP, CiP);
                    replace(newParent, newChild, iP);
                }

                iP.setHeight(newHeight);

                logq = Math.log(possibleDestinations);
            } else {
                iP.setHeight(newHeight);
                logq = 0.0;
            }
        }
        return logq;
    }

    private double getDelta() {
        if (!m_gaussian.get()) {
            return (Randomizer.nextDouble() * m_fSize) - (m_fSize / 2.0);
        } else {
            return Randomizer.nextGaussian() * m_fSize;
        }
    }

    private int intersectingEdges(Node node, double height, List<Node> directChildren) {
        final Node parent = node.getParent();

        if (parent.getHeight() < height) return 0;

        if (node.getHeight() < height) {
            if (directChildren != null) directChildren.add(node);
            return 1;
        }

        if (node.isLeaf()) {
            // TODO: verify that this makes sense
            return 0;
        } else {
            final int count = intersectingEdges(node.getLeft(), height, directChildren) +
                    intersectingEdges(node.getRight(), height, directChildren);
            return count;
        }
    }

    /**
     * automatic parameter tuning *
     */
    @Override
    public void optimize(final double logAlpha) {
        if (m_bOptimise.get()) {
            double fDelta = calcDelta(logAlpha);
            fDelta += Math.log(m_fSize);
            m_fSize = Math.exp(fDelta);
        }
    }

    @Override
    public double getCoercableParameterValue() {
        return m_fSize;
    }

    @Override
    public void setCoercableParameterValue(final double fValue) {
        m_fSize = fValue;
    }

    @Override
    public String getPerformanceSuggestion() {
        final double prob = m_nNrAccepted / (m_nNrAccepted + m_nNrRejected + 0.0);
        final double targetProb = getTargetAcceptanceProbability();

        double ratio = prob / targetProb;

        if (ratio > 2.0) ratio = 2.0;
        if (ratio < 0.5) ratio = 0.5;

        final double newDelta = m_fSize * ratio;

        final DecimalFormat formatter = new DecimalFormat("#.###");
        if (prob < 0.10) {
            return "Try decreasing size to about " + formatter.format(newDelta);
        } else if (prob > 0.40) {
            return "Try increasing size to about " + formatter.format(newDelta);
        } else return "";
    }

}
