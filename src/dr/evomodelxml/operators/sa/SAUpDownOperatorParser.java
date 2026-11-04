/*
 * SAUpDownOperatorParser.java
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

package dr.evomodelxml.operators.sa;

import dr.evomodel.operators.sa.*;
import dr.evomodel.tree.TreeModel;
import dr.evolution.util.TaxonList;
import dr.inference.model.Parameter;
import dr.inference.operators.MCMCOperator;
import dr.xml.*;

import java.util.ArrayList;
import java.util.List;

/** XML parser for SAUpDownOperator (port of BEAST 2's SA package operator; see the class). */
public class SAUpDownOperatorParser extends AbstractXMLObjectParser {

    public String getParserName() {
        return "saUpDownOperator";
    }

    public Class getReturnType() {
        return SAUpDownOperator.class;
    }

    public String getParserDescription() {
        return "Sampled-ancestor tree: scale a tree and parameters up and parameters down (BEAST 2 UpDownOperator)";
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        TreeModel tree = null;
        List<Parameter> up = new ArrayList<>(), down = new ArrayList<>();
        if (xo.hasChildNamed("up")) {
            XMLObject cxo = xo.getChild("up");
            for (int i = 0; i < cxo.getChildCount(); i++) {
                Object c = cxo.getChild(i);
                if (c instanceof TreeModel) {
                    tree = (TreeModel) c;
                } else {
                    up.add((Parameter) c);
                }
            }
        }
        if (xo.hasChildNamed("down")) {
            XMLObject cxo = xo.getChild("down");
            for (int i = 0; i < cxo.getChildCount(); i++) {
                down.add((Parameter) cxo.getChild(i));
            }
        }
        return new SAUpDownOperator(tree, up, down, xo.getDoubleAttribute("scaleFactor"),
                xo.getDoubleAttribute(MCMCOperator.WEIGHT));
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[]{AttributeRule.newDoubleRule(MCMCOperator.WEIGHT),
                AttributeRule.newDoubleRule("scaleFactor"),
                new ElementRule("up", new XMLSyntaxRule[]{new ElementRule(TreeModel.class, 0, 1),
                        new ElementRule(Parameter.class, 0, Integer.MAX_VALUE)}, true),
                new ElementRule("down", new XMLSyntaxRule[]{new ElementRule(Parameter.class, 1, Integer.MAX_VALUE)}, true)};
    }
}
