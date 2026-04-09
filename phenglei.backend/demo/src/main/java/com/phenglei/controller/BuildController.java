package com.phenglei.controller;

import com.phenglei.service.CmakeBuildService;
import lombok.RequiredArgsConstructor;
import org.springframework.web.bind.annotation.GetMapping;
import org.springframework.web.bind.annotation.RestController;

import java.util.HashMap;
import java.util.Map;

@RestController
@RequiredArgsConstructor
public class BuildController {

    private final CmakeBuildService cmakeBuildService;

    /**
     * 编译 EXE
     */
    @GetMapping("/build-exe")
    public Map<String, Object> buildExe() {
        Map<String, Object> result = new HashMap<>();
        try {
            Boolean success = cmakeBuildService.buildExe().get();

            if (success) {
                result.put("code", 200);
                result.put("msg", "编译成功");
                result.put("exePath", cmakeBuildService.getExePath());
                result.put("dllPath", cmakeBuildService.getTecioDllPath());
            } else {
                result.put("code", 500);
                result.put("msg", "编译失败");
            }
        } catch (Exception e) {
            result.put("code", 500);
            result.put("msg", e.getMessage());
        }
        return result;
    }

    /**
     * ✅ 运行 EXE
     */
    @GetMapping("/run-exe")
    public Map<String, Object> runExe() {
        Map<String, Object> result = new HashMap<>();
        try {
            Boolean success = cmakeBuildService.runExe().get();

            result.put("code", success ? 200 : 500);
            result.put("msg", success ? "运行成功" : "运行失败");

        } catch (Exception e) {
            result.put("code", 500);
            result.put("msg", e.getMessage());
        }
        return result;
    }
}