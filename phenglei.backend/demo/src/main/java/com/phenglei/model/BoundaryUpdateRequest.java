package com.phenglei.model;

import java.util.Map;

public class BoundaryUpdateRequest {

    private Integer nBoundaryConditions;
    private Map<String, Map<String, Object>> boundaries;

    public Integer getNBoundaryConditions() {
        return nBoundaryConditions;
    }

    public void setNBoundaryConditions(Integer nBoundaryConditions) {
        this.nBoundaryConditions = nBoundaryConditions;
    }

    public Map<String, Map<String, Object>> getBoundaries() {
        return boundaries;
    }

    public void setBoundaries(Map<String, Map<String, Object>> boundaries) {
        this.boundaries = boundaries;
    }
}