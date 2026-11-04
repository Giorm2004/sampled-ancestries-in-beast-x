/*
 * SampledAncestorOperatorTest.java
 *
 * Copyright © 2002-2024 the BEAST Development Team
 * http://beast.community/about
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
 *
 */

package test.dr.evomodel.operators;

import dr.app.beast.BeastParser;
import dr.evolution.tree.NodeRef;
import dr.evomodel.tree.SampledAncestorTrees;
import dr.evomodel.tree.TreeModel;
import dr.inference.operators.MCMCOperator;
import dr.inference.operators.OperatorSchedule;
import dr.math.MathUtils;
import dr.xml.XMLObject;
import junit.framework.TestCase;

import java.io.FileReader;
import java.util.Map;

/**
 * Structural invariants of the sampled-ancestor operators (examples/TestXML/testSampledAncestorOperators.xml): applying
 * every operator many times (each proposal kept, accepted or not, which is harsher than MCMC) always leaves a valid
 * tree -- no node above its parent, every sampled ancestor exactly at its fake parent's height, tips at the present
 * unmoved -- and the operators reach both trees with and without sampled ancestors. Distributional agreement with
 * BEAST 2 is documented in the XML.
 */
public class SampledAncestorOperatorTest extends TestCase {

    private static final String XML = "examples/TestXML/testSampledAncestorOperators.xml";

    public void testInvariants() throws Exception {
        MathUtils.setSeed(7);
        BeastParser parser = new BeastParser(new String[]{XML}, null, false, false, false, null);
        Map<String, XMLObject> store = parser.parse(new FileReader(XML), false);
        TreeModel tree = (TreeModel) store.get("treeModel").getNativeObject();
        OperatorSchedule schedule = (OperatorSchedule) store.get("operators").getNativeObject();
        int maxSA = 0, minSA = Integer.MAX_VALUE, valid = 0;
        for (int step = 0; step < 200000; step++) {
            MCMCOperator op = schedule.getOperator(schedule.getNextOperatorIndex());
            double logHR = op.operate();
            op.accept(0.0);                   // keep the proposal (operate() requires an accept or reject in between)
            if (logHR != Double.NEGATIVE_INFINITY) {
                valid++;
            }
            for (int k = 0; k < tree.getNodeCount(); k++) {
                NodeRef node = tree.getNode(k);
                if (!tree.isRoot(node)) {
                    assertTrue(op.getOperatorName() + ": node above its parent",
                            tree.getNodeHeight(node) <= tree.getNodeHeight(tree.getParent(node)));
                }
                if (tree.isExternal(node) && tree.getTaxonId(k).startsWith("m")) {
                    assertEquals(op.getOperatorName() + ": a tip at the present moved", 0.0, tree.getNodeHeight(node));
                }
            }
            int sa = SampledAncestorTrees.countSampledAncestors(tree);
            for (int k = 0; k < tree.getInternalNodeCount(); k++) {
                NodeRef node = tree.getInternalNode(k);
                if (SampledAncestorTrees.isFake(tree, node)) {
                    assertEquals(1, (SampledAncestorTrees.isSampledAncestor(tree, tree.getChild(node, 0)) ? 1 : 0)
                            + (SampledAncestorTrees.isSampledAncestor(tree, tree.getChild(node, 1)) ? 1 : 0));
                }
            }
            maxSA = Math.max(maxSA, sa);
            minSA = Math.min(minSA, sa);
        }
        System.out.println("200000 proposals, " + valid + " valid; sampled ancestors ranged " + minSA + ".." + maxSA);
        assertEquals(0, minSA);
        assertTrue(maxSA >= 3);
    }
}
