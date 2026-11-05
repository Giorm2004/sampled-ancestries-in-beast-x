/*
 * TimeDependentBranchRates.java
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

package dr.evomodel.branchratemodel;

import dr.evolution.tree.BranchRates;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;

/**
 * Branch rates that vary along a branch with calendar time. A branch's length in substitutions is then
 * b = integral of r(t) over the branch, not (average rate) x (duration), so
 *
 *     db / d h_parent = r(h_parent),     db / d h_child = -r(h_child),
 *
 * where r(t) is the instantaneous rate of that branch at time t. The node-height gradient of a tree likelihood
 * must use these, not getBranchRate (the average), whenever the rate depends on time.
 *
 * @author Giovanni Mazzeo
 */
public interface TimeDependentBranchRates extends BranchRates {

    /** The instantaneous rate on the branch above node at a time between the heights of node and its parent. */
    double getBranchRateAtTime(Tree tree, NodeRef node, double time);
}
