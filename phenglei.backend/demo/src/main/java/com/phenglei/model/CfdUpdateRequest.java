package com.phenglei.model;

public class CfdUpdateRequest {

    private Integer maxSimuStep;
    private Integer intervalStepFlow;
    private Integer intervalStepPlot;
    private Integer intervalStepForce;
    private Integer intervalStepRes;
    private Integer ifLowSpeedPrecon;

    private Double refMachNumber;
    private Double attackd;
    private Double angleSlide;
    private Integer inflowParaType;
    private Integer refReNumber;
    private Double refDimensionalTemperature;
    private Double gridScaleFactor;

    private Double forceReferenceLengthSpanWise;
    private Double forceReferenceLength;
    private Double forceReferenceArea;
    private Double TorqueRefX;
    private Double TorqueRefY;
    private Double TorqueRefZ;

    private Integer viscousType;
    private String viscousName;
    private Integer DESType;

    private Integer roeEntropyFixMethod;
    private Double roeEntropyScale;

    private String str_limiter_name;
    private String uns_limiter_name;
    private Double venkatCoeff;

    private Integer iunsteady;
    private Double CFLEnd;
    private Integer nLUSGSSweeps;
    private Integer nMGLevel;
    private Integer flowInitStep;

    private String gridfile;
    private Integer plotFieldType;
    private Integer nVisualVariables;
    private String visualVariables;

    private Integer reconmeth;
    private Integer limitVariables;
    private Integer limitVector;

    public Integer getMaxSimuStep() {
        return maxSimuStep;
    }

    public void setMaxSimuStep(Integer maxSimuStep) {
        this.maxSimuStep = maxSimuStep;
    }

    public Integer getIntervalStepFlow() {
        return intervalStepFlow;
    }

    public void setIntervalStepFlow(Integer intervalStepFlow) {
        this.intervalStepFlow = intervalStepFlow;
    }

    public Integer getIntervalStepPlot() {
        return intervalStepPlot;
    }

    public void setIntervalStepPlot(Integer intervalStepPlot) {
        this.intervalStepPlot = intervalStepPlot;
    }

    public Integer getIntervalStepForce() {
        return intervalStepForce;
    }

    public void setIntervalStepForce(Integer intervalStepForce) {
        this.intervalStepForce = intervalStepForce;
    }

    public Integer getIntervalStepRes() {
        return intervalStepRes;
    }

    public void setIntervalStepRes(Integer intervalStepRes) {
        this.intervalStepRes = intervalStepRes;
    }

    public Integer getIfLowSpeedPrecon() {
        return ifLowSpeedPrecon;
    }

    public void setIfLowSpeedPrecon(Integer ifLowSpeedPrecon) {
        this.ifLowSpeedPrecon = ifLowSpeedPrecon;
    }

    public Double getRefMachNumber() {
        return refMachNumber;
    }

    public void setRefMachNumber(Double refMachNumber) {
        this.refMachNumber = refMachNumber;
    }

    public Double getAttackd() {
        return attackd;
    }

    public void setAttackd(Double attackd) {
        this.attackd = attackd;
    }

    public Double getAngleSlide() {
        return angleSlide;
    }

    public void setAngleSlide(Double angleSlide) {
        this.angleSlide = angleSlide;
    }

    public Integer getInflowParaType() {
        return inflowParaType;
    }

    public void setInflowParaType(Integer inflowParaType) {
        this.inflowParaType = inflowParaType;
    }

    public Integer getRefReNumber() {
        return refReNumber;
    }

    public void setRefReNumber(Integer refReNumber) {
        this.refReNumber = refReNumber;
    }

    public Double getRefDimensionalTemperature() {
        return refDimensionalTemperature;
    }

    public void setRefDimensionalTemperature(Double refDimensionalTemperature) {
        this.refDimensionalTemperature = refDimensionalTemperature;
    }

    public Double getGridScaleFactor() {
        return gridScaleFactor;
    }

    public void setGridScaleFactor(Double gridScaleFactor) {
        this.gridScaleFactor = gridScaleFactor;
    }

    public Double getForceReferenceLengthSpanWise() {
        return forceReferenceLengthSpanWise;
    }

    public void setForceReferenceLengthSpanWise(Double forceReferenceLengthSpanWise) {
        this.forceReferenceLengthSpanWise = forceReferenceLengthSpanWise;
    }

    public Double getForceReferenceLength() {
        return forceReferenceLength;
    }

    public void setForceReferenceLength(Double forceReferenceLength) {
        this.forceReferenceLength = forceReferenceLength;
    }

    public Double getForceReferenceArea() {
        return forceReferenceArea;
    }

    public void setForceReferenceArea(Double forceReferenceArea) {
        this.forceReferenceArea = forceReferenceArea;
    }

    public Double getTorqueRefX() {
        return TorqueRefX;
    }

    public void setTorqueRefX(Double torqueRefX) {
        TorqueRefX = torqueRefX;
    }

    public Double getTorqueRefY() {
        return TorqueRefY;
    }

    public void setTorqueRefY(Double torqueRefY) {
        TorqueRefY = torqueRefY;
    }

    public Double getTorqueRefZ() {
        return TorqueRefZ;
    }

    public void setTorqueRefZ(Double torqueRefZ) {
        TorqueRefZ = torqueRefZ;
    }

    public Integer getViscousType() {
        return viscousType;
    }

    public void setViscousType(Integer viscousType) {
        this.viscousType = viscousType;
    }

    public String getViscousName() {
        return viscousName;
    }

    public void setViscousName(String viscousName) {
        this.viscousName = viscousName;
    }

    public Integer getDESType() {
        return DESType;
    }

    public void setDESType(Integer DESType) {
        this.DESType = DESType;
    }

    public Integer getRoeEntropyFixMethod() {
        return roeEntropyFixMethod;
    }

    public void setRoeEntropyFixMethod(Integer roeEntropyFixMethod) {
        this.roeEntropyFixMethod = roeEntropyFixMethod;
    }

    public Double getRoeEntropyScale() {
        return roeEntropyScale;
    }

    public void setRoeEntropyScale(Double roeEntropyScale) {
        this.roeEntropyScale = roeEntropyScale;
    }

    public String getStr_limiter_name() {
        return str_limiter_name;
    }

    public void setStr_limiter_name(String str_limiter_name) {
        this.str_limiter_name = str_limiter_name;
    }

    public String getUns_limiter_name() {
        return uns_limiter_name;
    }

    public void setUns_limiter_name(String uns_limiter_name) {
        this.uns_limiter_name = uns_limiter_name;
    }

    public Double getVenkatCoeff() {
        return venkatCoeff;
    }

    public void setVenkatCoeff(Double venkatCoeff) {
        this.venkatCoeff = venkatCoeff;
    }

    public Integer getIunsteady() {
        return iunsteady;
    }

    public void setIunsteady(Integer iunsteady) {
        this.iunsteady = iunsteady;
    }

    public Double getCFLEnd() {
        return CFLEnd;
    }

    public void setCFLEnd(Double CFLEnd) {
        this.CFLEnd = CFLEnd;
    }

    public Integer getNLUSGSSweeps() {
        return nLUSGSSweeps;
    }

    public void setNLUSGSSweeps(Integer nLUSGSSweeps) {
        this.nLUSGSSweeps = nLUSGSSweeps;
    }

    public Integer getNMGLevel() {
        return nMGLevel;
    }

    public void setNMGLevel(Integer NMGLevel) {
        this.nMGLevel = NMGLevel;
    }

    public Integer getFlowInitStep() {
        return flowInitStep;
    }

    public void setFlowInitStep(Integer flowInitStep) {
        this.flowInitStep = flowInitStep;
    }

    public String getGridfile() {
        return gridfile;
    }

    public void setGridfile(String gridfile) {
        this.gridfile = gridfile;
    }

    public Integer getPlotFieldType() {
        return plotFieldType;
    }

    public void setPlotFieldType(Integer plotFieldType) {
        this.plotFieldType = plotFieldType;
    }

    public Integer getNVisualVariables() {
        return nVisualVariables;
    }

    public void setNVisualVariables(Integer nVisualVariables) {
        this.nVisualVariables = nVisualVariables;
    }

    public String getVisualVariables() {
        return visualVariables;
    }

    public void setVisualVariables(String visualVariables) {
        this.visualVariables = visualVariables;
    }

    public Integer getReconmeth() {
        return reconmeth;
    }

    public void setReconmeth(Integer reconmeth) {
        this.reconmeth = reconmeth;
    }

    public Integer getLimitVariables() {
        return limitVariables;
    }

    public void setLimitVariables(Integer limitVariables) {
        this.limitVariables = limitVariables;
    }

    public Integer getLimitVector() {
        return limitVector;
    }

    public void setLimitVector(Integer limitVector) {
        this.limitVector = limitVector;
    }
}