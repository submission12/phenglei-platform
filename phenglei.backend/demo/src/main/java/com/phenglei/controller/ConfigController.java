package com.phenglei.controller;

import com.phenglei.model.BoundaryUpdateRequest;
import com.phenglei.model.CfdUpdateRequest;
import com.phenglei.model.KeyUpdateRequest;
import com.phenglei.service.HyparaConfigService;
import com.phenglei.service.ConfigResetService;

import org.springframework.http.HttpStatus;
import org.springframework.web.bind.annotation.*;
import org.springframework.web.server.ResponseStatusException;

import java.util.Map;

@RestController
@RequestMapping("/config")
public class ConfigController {

    private final HyparaConfigService hyparaConfigService;
    private final ConfigResetService resetService; // ⭐ 新增

    public ConfigController(HyparaConfigService hyparaConfigService,
                            ConfigResetService resetService) {
        this.hyparaConfigService = hyparaConfigService;
        this.resetService = resetService;
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

    // ✅ ⭐ 新增：恢复默认配置
    @PostMapping("/reset")
    public Map<String, Object> reset() {
        try {
            resetService.resetToDefault();
            return Map.of(
                    "ok", true,
                    "msg", "已恢复默认配置"
            );
        } catch (Exception e) {
            throw new ResponseStatusException(HttpStatus.INTERNAL_SERVER_ERROR, "恢复失败", e);
        }
    }
}