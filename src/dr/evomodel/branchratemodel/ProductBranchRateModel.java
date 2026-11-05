/*
 * ProductBranchRateModel.java
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

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.xml.*;

import java.util.ArrayList;
import java.util.List;

/**
 * Differentiable product of branch-rate models, r_b = prod_k F_k,b: e.g. a relaxed clock (per-branch rates,
 * ArbitraryBranchRates) times a time-varying clock (TimeVaryingBranchRateModel). At most one factor may depend
 * on time; the branch length in substitutions is then (prod of the other factors) x integral of that factor's
 * rate, which is what makes the product exact and the node-height chain rule simple:
 *
 *     db/dh = (prod_{j != td} F_j,b) x r_td(h)    (TimeDependentBranchRates).
 *
 * Gradients w.r.t. a factor's parameters go through a Factor view (XML productBranchRateFactor), which a
 * branchRateGradient takes alongside the tree likelihood: dlogL/dF_k,b = dlogL/dr_b x prod_{j != k} F_j,b,
 * then the factor's own updateGradientLogDensity.
 *
 * @author Giovanni Mazzeo
 */
public class ProductBranchRateModel extends AbstractBranchRateModel
        implements DifferentiableBranchRates, TimeDependentBranchRates {

    private final Tree tree;
    private final List<DifferentiableBranchRates> factors;
    private final int timeDependent;      // index of the time-dependent factor, -1 if none

    public ProductBranchRateModel(Tree tree, List<DifferentiableBranchRates> factors) {
        super("productBranchRates");
        this.tree = tree;
        this.factors = new ArrayList<>(factors);
        int td = -1;
        for (int k = 0; k < factors.size(); k++) {
            DifferentiableBranchRates f = factors.get(k);
            if (!(f instanceof Model)) {
                throw new IllegalArgumentException("Factor " + k + " is not a model");
            }
            addModel((Model) f);
            if (f instanceof TimeDependentBranchRates) {
                if (td >= 0) {
                    throw new IllegalArgumentException("At most one factor may depend on time (the product of two "
                            + "branch averages is not the average of the product)");
                }
                td = k;
            }
        }
        this.timeDependent = td;
    }

    public List<DifferentiableBranchRates> getFactors() {
        return factors;
    }

    private double productExcept(Tree tree, NodeRef node, int skip) {
        double r = 1.0;
        for (int k = 0; k < factors.size(); k++) {
            if (k != skip) {
                r *= factors.get(k).getBranchRate(tree, node);
            }
        }
        return r;
    }

    @Override
    public double getBranchRate(Tree tree, NodeRef node) {
        return productExcept(tree, node, -1);
    }

    @Override
    public double getBranchRateAtTime(Tree tree, NodeRef node, double time) {
        if (timeDependent < 0) {
            return getBranchRate(tree, node);
        }
        return productExcept(tree, node, timeDependent)
                * ((TimeDependentBranchRates) factors.get(timeDependent)).getBranchRateAtTime(tree, node, time);
    }

    @Override
    public Tree getTree() {
        return tree;
    }

    @Override
    public int getParameterIndexFromNode(NodeRef node) {
        int n = node.getNumber();
        return n > tree.getRoot().getNumber() ? n - 1 : n;
    }

    private static UnsupportedOperationException useFactor() {
        return new UnsupportedOperationException("Take gradients w.r.t. a factor: use a productBranchRateFactor");
    }

    @Override
    public double getBranchRateDifferential(Tree tree, NodeRef node) {
        throw useFactor();
    }

    @Override
    public double getBranchRateSecondDifferential(Tree tree, NodeRef node) {
        throw useFactor();
    }

    @Override
    public Parameter getRateParameter() {
        throw useFactor();
    }

    @Override
    public ArbitraryBranchRates.BranchRateTransform getTransform() {
        throw useFactor();
    }

    @Override
    public double[] updateGradientLogDensity(double[] gradient, double[] value, int from, int to) {
        throw useFactor();
    }

    @Override
    public double[] updateDiagonalHessianLogDensity(double[] diagonalHessian, double[] gradient, double[] value,
                                                    int from, int to) {
        throw useFactor();
    }

    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {
        fireModelChanged();
    }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) { }

    @Override
    protected void storeState() { }

    @Override
    protected void restoreState() { }

    @Override
    protected void acceptState() { }

    /** Factor k seen through the product: its parameters, with the chain rule through the other factors. */
    public class Factor implements DifferentiableBranchRates {

        private final int k;
        private final DifferentiableBranchRates factor;

        public Factor(int k) {
            this.k = k;
            this.factor = factors.get(k);
        }

        public ProductBranchRateModel getProduct() {
            return ProductBranchRateModel.this;
        }

        @Override
        public Tree getTree() {
            return tree;
        }

        @Override
        public double getBranchRate(Tree tree, NodeRef node) {
            return factor.getBranchRate(tree, node);
        }

        @Override
        public double getBranchRateDifferential(Tree tree, NodeRef node) {
            return factor.getBranchRateDifferential(tree, node) * productExcept(tree, node, k);
        }

        @Override
        public double getBranchRateSecondDifferential(Tree tree, NodeRef node) {
            return factor.getBranchRateSecondDifferential(tree, node) * productExcept(tree, node, k);
        }

        @Override
        public Parameter getRateParameter() {
            return factor.getRateParameter();
        }

        @Override
        public int getParameterIndexFromNode(NodeRef node) {
            return factor.getParameterIndexFromNode(node);
        }

        @Override
        public ArbitraryBranchRates.BranchRateTransform getTransform() {
            return factor.getTransform();
        }

        @Override
        public double[] updateGradientLogDensity(double[] gradient, double[] value, int from, int to) {
            return factor.updateGradientLogDensity(gradient, value, from, to);
        }

        @Override
        public double[] updateDiagonalHessianLogDensity(double[] diagonalHessian, double[] gradient, double[] value,
                                                        int from, int to) {
            return factor.updateDiagonalHessianLogDensity(diagonalHessian, gradient, value, from, to);
        }
    }

    // **************************************************************
    // XML
    // **************************************************************

    public static final XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return "productBranchRates";
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {
            Tree tree = (Tree) xo.getChild(Tree.class);
            List<DifferentiableBranchRates> factors = new ArrayList<>();
            for (int i = 0; i < xo.getChildCount(); i++) {
                Object child = xo.getChild(i);
                if (child instanceof DifferentiableBranchRates) {
                    factors.add((DifferentiableBranchRates) child);
                }
            }
            if (factors.size() < 2) {
                throw new XMLParseException("productBranchRates needs at least two differentiable branch-rate models");
            }
            try {
                return new ProductBranchRateModel(tree, factors);
            } catch (IllegalArgumentException e) {
                throw new XMLParseException(e.getMessage());
            }
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return new XMLSyntaxRule[]{
                    new ElementRule(Tree.class),
                    new ElementRule(DifferentiableBranchRates.class, 2, Integer.MAX_VALUE)
            };
        }

        public String getParserDescription() {
            return "Differentiable product of branch-rate models (at most one time-dependent), e.g. relaxed x time-varying";
        }

        public Class getReturnType() {
            return ProductBranchRateModel.class;
        }
    };

    public static final XMLObjectParser FACTOR_PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return "productBranchRateFactor";
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {
            ProductBranchRateModel product = (ProductBranchRateModel) xo.getChild(ProductBranchRateModel.class);
            for (int i = 0; i < xo.getChildCount(); i++) {
                Object child = xo.getChild(i);
                if (child != product) {
                    int k = product.factors.indexOf(child);
                    if (k < 0) {
                        throw new XMLParseException("Not a factor of " + product.getId() + ": " + child);
                    }
                    return product.new Factor(k);
                }
            }
            throw new XMLParseException("Name the factor: productBranchRateFactor needs the product and one of its factors");
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return new XMLSyntaxRule[]{
                    new ElementRule(ProductBranchRateModel.class),
                    new ElementRule(DifferentiableBranchRates.class, 1, 2)
            };
        }

        public String getParserDescription() {
            return "One factor of a productBranchRates, for branchRateGradient";
        }

        public Class getReturnType() {
            return Factor.class;
        }
    };
}
