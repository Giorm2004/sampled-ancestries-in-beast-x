/*
 * NewBirthDeathSerialSamplingModelParser.java
 *
 * Copyright © 2002-2026 the BEAST Development Team
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

package dr.evomodelxml.birthdeath;

import dr.evolution.util.Units;
import dr.evomodel.birthdeath.EpisodicBirthDeathSamplingModel;
import dr.evomodel.birthdeath.RemovalAwareBirthDeathRates;
import dr.evoxml.util.XMLUnits;
import dr.inference.model.Parameter;
import dr.xml.*;

import java.util.Arrays;
import java.util.logging.Logger;

/**
 * Parser for the EpisodicBirthDeathSamplingModel (EBDS).
 */
public class EpisodicBirthDeathSamplingModelParser extends AbstractXMLObjectParser {

    public static final String EPISODIC_BIRTH_DEATH_SAMPLING_MODEL = "episodicBirthDeathSamplingModel";
    public static final String BIRTH_RATE = "birthRate";
    public static final String DEATH_RATE = "deathRate";
    public static final String SAMPLING_RATE = "samplingRate";
    public static final String TREATMENT_PROBABILITY = "treatmentProbability";
    public static final String SAMPLING_PROBABILITY = "samplingProbability";
    public static final String ORIGIN = "origin";
    public static final String CONDITION = "conditionOnSurvival";
    public static final String NUM_GRID_POINTS = "numGridPoints";
    public static final String CUT_OFF = "cutOff";
    public static final String R0 = "R0";
    public static final String D = "D";
    public static final String S = "S";
    public static final String GRADIENT_FLAG = "gradientFlag";
    public static final String GRIDS = "grids";
    public static final String BREAK_POINTS = "breakPoints";
    public static final String REPRODUCTIVE_NUMBER = "reproductiveNumber";
    public static final String BECOME_UNINFECTIOUS_RATE = "becomeUninfectiousRate";
    public static final String SAMPLING_PROPORTION = "samplingProportion";
    public static final String SAMPLED_ANCESTORS = "sampledAncestors";

    public String getParserName() {
        return EPISODIC_BIRTH_DEATH_SAMPLING_MODEL;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        Units.Type units = XMLUnits.Utils.getUnitsAttr(xo);
        final Parameter originParameter = (Parameter) xo.getChild(ORIGIN).getChild(Parameter.class);
        final Parameter numGridPoints = xo.hasChildNamed(NUM_GRID_POINTS) ? (Parameter) xo.getElementFirstChild(NUM_GRID_POINTS): new Parameter.Default(1.0);
        final Parameter cutoff = xo.hasChildNamed(CUT_OFF) ? (Parameter) xo.getElementFirstChild(CUT_OFF) : new Parameter.Default(Double.POSITIVE_INFINITY);
        final Parameter grids = xo.hasChildNamed(GRIDS) ? (Parameter) xo.getChild(GRIDS).getChild(Parameter.class) : null;
        final Parameter breakPoints = xo.hasChildNamed(BREAK_POINTS) ? (Parameter) xo.getChild(BREAK_POINTS).getChild(Parameter.class) : null;
        if (breakPoints != null && grids != null) {
            throw new XMLParseException("give either " + GRIDS + " or " + BREAK_POINTS + ", not both");
        }
        final int numIntervals = breakPoints != null ? breakPoints.getDimension() + 1 : (int) numGridPoints.getParameterValue(0);
        if (breakPoints != null && xo.hasChildNamed(NUM_GRID_POINTS) && (int) numGridPoints.getParameterValue(0) != numIntervals) {
            throw new XMLParseException(NUM_GRID_POINTS + " must be the dimension of " + BREAK_POINTS + " plus one");
        }

        Parameter lambda = null;
        Parameter mu = null;
        Parameter psi = null;
        Parameter r = null;
        Parameter rho = null;

        Parameter R0Parameter = null;
        Parameter DParameter = null;
        Parameter SParameter = null;

        boolean[] gradientFlags = new boolean[5];
        Arrays.fill(gradientFlags, true);

        RemovalAwareBirthDeathRates epiRates = null;

        if (xo.hasChildNamed(REPRODUCTIVE_NUMBER)) {
            // BDSKY's epidemiological parameterisation with removal probability (see RemovalAwareBirthDeathRates)
            epiRates = new RemovalAwareBirthDeathRates(
                    (Parameter) xo.getChild(REPRODUCTIVE_NUMBER).getChild(Parameter.class),
                    (Parameter) xo.getChild(BECOME_UNINFECTIOUS_RATE).getChild(Parameter.class),
                    (Parameter) xo.getChild(SAMPLING_PROPORTION).getChild(Parameter.class),
                    (Parameter) xo.getChild(TREATMENT_PROBABILITY).getChild(Parameter.class),
                    numIntervals);
            lambda = epiRates.birthRate;
            mu = epiRates.deathRate;
            psi = epiRates.samplingRate;
            Arrays.fill(gradientFlags, 0, 3, false);       // no gradients wrt R0, delta, s

        } else if (xo.hasChildNamed(BIRTH_RATE)) {

            lambda = (Parameter) xo.getChild(BIRTH_RATE).getChild(Parameter.class);
            mu = (Parameter) xo.getChild(DEATH_RATE).getChild(Parameter.class);
            psi = (Parameter) xo.getChild(SAMPLING_RATE).getChild(Parameter.class);

            if (xo.getChild(BIRTH_RATE).hasAttribute(GRADIENT_FLAG)) {
                gradientFlags[0] = xo.getChild(BIRTH_RATE).getAttribute(GRADIENT_FLAG, true);
            }
            if (xo.getChild(DEATH_RATE).hasAttribute(GRADIENT_FLAG)) {
                gradientFlags[1] = xo.getChild(DEATH_RATE).getAttribute(GRADIENT_FLAG, true);
            }
            if (xo.getChild(SAMPLING_RATE).hasAttribute(GRADIENT_FLAG)) {
                gradientFlags[2] = xo.getChild(SAMPLING_RATE).getAttribute(GRADIENT_FLAG, true);
            }
        } else {

            R0Parameter = (Parameter) xo.getChild(R0).getChild(Parameter.class);
            DParameter = (Parameter) xo.getChild(D).getChild(Parameter.class);
            SParameter = (Parameter) xo.getChild(S).getChild(Parameter.class);


            if (xo.getChild(R0).hasAttribute(GRADIENT_FLAG)) {
                gradientFlags[0] = xo.getChild(R0).getAttribute(GRADIENT_FLAG, true);
            }
            if (xo.getChild(D).hasAttribute(GRADIENT_FLAG)) {
                gradientFlags[1] = xo.getChild(D).getAttribute(GRADIENT_FLAG, true);
            }
            if (xo.getChild(S).hasAttribute(GRADIENT_FLAG)) {
                gradientFlags[2] = xo.getChild(S).getAttribute(GRADIENT_FLAG, true);
        }
        }

        r = (Parameter) xo.getChild(TREATMENT_PROBABILITY).getChild(Parameter.class);
        rho = (Parameter) xo.getChild(SAMPLING_PROBABILITY).getChild(Parameter.class);

        if (xo.getChild(TREATMENT_PROBABILITY).hasAttribute(GRADIENT_FLAG)) {
            gradientFlags[3] = xo.getChild(TREATMENT_PROBABILITY).getAttribute(GRADIENT_FLAG, true);
        }
        if (xo.getChild(SAMPLING_PROBABILITY).hasAttribute(GRADIENT_FLAG)) {
            gradientFlags[4] = xo.getChild(SAMPLING_PROBABILITY).getAttribute(GRADIENT_FLAG, true);
        }

        boolean condition = xo.getAttribute(CONDITION, false);

        EpisodicBirthDeathSamplingModel model;

        if (R0Parameter != null) {
            model = EpisodicBirthDeathSamplingModel.createWithCompoundParameters(
                    R0Parameter, DParameter, SParameter,
                    r, rho, originParameter,
                    condition, numIntervals,
                    cutoff.getParameterValue(0), units);
        } else {
            model = new EpisodicBirthDeathSamplingModel(lambda, mu, psi, r, rho, originParameter,
                    condition, numIntervals,
                    cutoff.getParameterValue(0), units);
        }


        model.setupGradientFlags(gradientFlags);
        model.setSampledAncestors(xo.getAttribute(SAMPLED_ANCESTORS, false));

        if (breakPoints != null) {
            model.setupTimeline(breakPoints);
        } else {
            model.setupTimeline(grids != null ? grids.getParameterValues() : null);
        }

        Logger.getLogger("dr.evomodel").info("Using birth-death serial sampling model");
        Logger.getLogger("dr.evomodel").info("\tCondition on root: " + condition);
        Logger.getLogger("dr.evomodel").info("\tGrid size: " + numIntervals);
        Logger.getLogger("dr.evomodel").info("\tGrid end: " + cutoff.getParameterValue(0));
        Logger.getLogger("dr.evomodel").info("\tParameter mode: " + (epiRates != null ? "epidemiological (BDSKY, removal-aware)"
                : R0Parameter != null ? "compound" : "raw"));
        if (breakPoints != null) {
            Logger.getLogger("dr.evomodel").info("\tBreak points: parameter " + breakPoints.getId());
        }
        if (model.hasSampledAncestors()) {
            Logger.getLogger("dr.evomodel").info("\tSampled ancestors allowed (BDSKY terms, labelling log 2 per non-fake internal node)");
        }

        return model;
    }

    public String getParserDescription() {
        return "A birth-death model with serial sampling through time and sampled ancestors.";
    }

    public Class getReturnType() {
        return EpisodicBirthDeathSamplingModel.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newBooleanRule(CONDITION, true),
            AttributeRule.newBooleanRule(SAMPLED_ANCESTORS, true),
            new XORRule(new XMLSyntaxRule[]{
                    new ElementRule(BIRTH_RATE,
                            new XMLSyntaxRule[]{
                                new ElementRule(Parameter.class),
                                AttributeRule.newBooleanRule(GRADIENT_FLAG, true)
                            }),
                    new ElementRule(R0,
                            new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
                    new ElementRule(REPRODUCTIVE_NUMBER,
                            new XMLSyntaxRule[]{new ElementRule(Parameter.class)})
            }),
            new ElementRule(BECOME_UNINFECTIOUS_RATE,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}, true),
            new ElementRule(SAMPLING_PROPORTION,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}, true),
            new ElementRule(BREAK_POINTS,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}, true),
            new ElementRule(DEATH_RATE,
                    new XMLSyntaxRule[]{
                        new ElementRule(Parameter.class),
                        AttributeRule.newBooleanRule(GRADIENT_FLAG, true)
                    }, true),
            new ElementRule(SAMPLING_RATE,
                    new XMLSyntaxRule[]{
                        new ElementRule(Parameter.class),
                        AttributeRule.newBooleanRule(GRADIENT_FLAG, true)
                    }, true),
            new ElementRule(D,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}, true),
            new ElementRule(S,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}, true),
            new ElementRule(TREATMENT_PROBABILITY,
                    new XMLSyntaxRule[]{
                        new ElementRule(Parameter.class),
                        AttributeRule.newBooleanRule(GRADIENT_FLAG, true)
                    }),
            new ElementRule(SAMPLING_PROBABILITY,
                    new XMLSyntaxRule[]{
                        new ElementRule(Parameter.class),
                        AttributeRule.newBooleanRule(GRADIENT_FLAG, true)
                    }),
            new ElementRule(ORIGIN,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
            new ElementRule(NUM_GRID_POINTS,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)},true),
            new ElementRule(CUT_OFF,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)},true),
            new ElementRule(GRIDS,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}, true),
            XMLUnits.SYNTAX_RULES[0]
    };
}