package com.phenglei.model;

public class KeyUpdateRequest {
    private Integer ndim;
    private Integer nparafile;
    private Integer nsimutask;
    private String parafilename;

    public Integer getNdim() {
        return ndim;
    }

    public void setNdim(Integer ndim) {
        this.ndim = ndim;
    }

    public Integer getNparafile() {
        return nparafile;
    }

    public void setNparafile(Integer nparafile) {
        this.nparafile = nparafile;
    }

    public Integer getNsimutask() {
        return nsimutask;
    }

    public void setNsimutask(Integer nsimutask) {
        this.nsimutask = nsimutask;
    }

    public String getParafilename() {
        return parafilename;
    }

    public void setParafilename(String parafilename) {
        this.parafilename = parafilename;
    }
}