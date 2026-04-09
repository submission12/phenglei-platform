package com.phenglei.controller;

import com.phenglei.model.BoundaryUpdateRequest;
import com.phenglei.model.CfdUpdateRequest;
import com.phenglei.model.KeyUpdateRequest;
import com.phenglei.service.HyparaConfigService;
import org.springframework.http.HttpStatus;
import org.springframework.web.bind.annotation.GetMapping;
import org.springframework.web.bind.annotation.PatchMapping;
import org.springframework.web.bind.annotation.RequestBody;
import org.springframework.web.bind.annotation.RequestMapping;
import org.springframework.web.bind.annotation.RestController;
import org.springframework.web.server.ResponseStatusException;

import java.util.Map;

@RestController
@RequestMapping("/config")
public class ConfigController {

    private final HyparaConfigService hyparaConfigService;

    public ConfigController(HyparaConfigService hyparaConfigService) {
        this.hyparaConfigService = hyparaConfigService;
    }

    @GetMapping("/key")
    public Map<String, Object> getKeyConfig() {
        try {
            return hyparaConfigService.getKeyConfig();
        } catch (Exception e) {
            throw new ResponseStatusException(HttpStatus.INTERNAL_SERVER_ERROR, e.getMessage(), e);
        }
    }

    @PatchMapping("/key")
    public Map<String, Object> updateKeyConfig(@RequestBody KeyUpdateRequest request) {
        try {
            return hyparaConfigService.updateKeyConfig(request);
        } catch (IllegalArgumentException e) {
            throw new ResponseStatusException(HttpStatus.BAD_REQUEST, e.getMessage(), e);
        } catch (Exception e) {
            throw new ResponseStatusException(HttpStatus.INTERNAL_SERVER_ERROR, e.getMessage(), e);
        }
    }

    @GetMapping("/cfd")
    public Map<String, Object> getCfdConfig() {
        try {
            return hyparaConfigService.getCfdConfig();
        } catch (Exception e) {
            throw new ResponseStatusException(HttpStatus.INTERNAL_SERVER_ERROR, e.getMessage(), e);
        }
    }

    @PatchMapping("/cfd")
    public Map<String, Object> updateCfdConfig(@RequestBody CfdUpdateRequest request) {
        try {
            return hyparaConfigService.updateCfdConfig(request);
        } catch (IllegalArgumentException e) {
            throw new ResponseStatusException(HttpStatus.BAD_REQUEST, e.getMessage(), e);
        } catch (Exception e) {
            throw new ResponseStatusException(HttpStatus.INTERNAL_SERVER_ERROR, e.getMessage(), e);
        }
    }

    @GetMapping("/boundary")
    public Map<String, Object> getBoundaryConfig() {
        try {
            return hyparaConfigService.getBoundaryConfig();
        } catch (Exception e) {
            throw new ResponseStatusException(HttpStatus.INTERNAL_SERVER_ERROR, e.getMessage(), e);
        }
    }

    @PatchMapping("/boundary")
    public Map<String, Object> updateBoundaryConfig(@RequestBody BoundaryUpdateRequest request) {
        try {
            return hyparaConfigService.updateBoundaryConfig(request);
        } catch (IllegalArgumentException e) {
            throw new ResponseStatusException(HttpStatus.BAD_REQUEST, e.getMessage(), e);
        } catch (Exception e) {
            throw new ResponseStatusException(HttpStatus.INTERNAL_SERVER_ERROR, e.getMessage(), e);
        }
    }
}