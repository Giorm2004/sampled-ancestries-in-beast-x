/*
 * SAScaleOperatorParser.java
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

/** XML parser for SAScaleOperator (port of BEAST 2's SA package operator; see the class). */
public class SAScaleOperatorParser extends AbstractXMLObjectParser {

    public String getParserName() {
        return "saScaleOperator";
    }

    public Class getReturnType() {
        return SAScaleOperator.class;
    }

    public String getParserDescription() {
        return "Sampled-ancestor tree: scale non-fake internal node heights, or the root only (BEAST 2 SA SAScaleOperator)";
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        return new SAScaleOperator((TreeModel) xo.getChild(TreeModel.class), xo.getDoubleAttribute("scaleFactor"),
                xo.getAttribute("rootOnly", false), xo.getDoubleAttribute(MCMCOperator.WEIGHT));
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[]{AttributeRule.newDoubleRule(MCMCOperator.WEIGHT),
                AttributeRule.newDoubleRule("scaleFactor"), AttributeRule.newBooleanRule("rootOnly", true),
                new ElementRule(TreeModel.class)};
    }
}
