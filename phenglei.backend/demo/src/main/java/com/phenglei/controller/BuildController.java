package com.phenglei.controller;

import com.phenglei.service.CmakeBuildService;
import com.phenglei.service.LogService;
import lombok.RequiredArgsConstructor;
import org.springframework.web.bind.annotation.GetMapping;
import org.springframework.web.bind.annotation.RestController;

import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.util.HashMap;
import java.util.Map;

@RestController
@RequiredArgsConstructor
public class BuildController {

    private final CmakeBuildService cmakeBuildService;
    private final LogService logService;

    /**
     * 编译
     */
    @GetMapping("/build-exe")
    public Map<String, Object> buildExe() {

        Map<String, Object> result = new HashMap<>();

        try {
            boolean ok = cmakeBuildService.buildExe().get();

            result.put("code", ok ? 200 : 500);
            result.put("msg", ok ? "编译成功" : "编译失败");

        } catch (Exception e) {
            result.put("code", 500);
            result.put("msg", e.getMessage());
        }

        return result;
    }

    /**
     * ⭐ 运行 EXE + SSE实时日志
     */
    @GetMapping("/run-exe")
    public Map<String, Object> runExe() {

        Map<String, Object> result = new HashMap<>();
        String taskId = "run";

        try {
            Process process = cmakeBuildService.startExe();

            new Thread(() -> {

                try (BufferedReader reader =
                             new BufferedReader(
                                     new InputStreamReader(process.getInputStream())
                             )) {

                    String line;

                    while ((line = reader.readLine()) != null) {
                        logService.sendLog(taskId, line);
                    }

                    int code = process.waitFor();

                    logService.sendLog(taskId, "结束 exitCode=" + code);
                    logService.sendProgress(taskId, 100);
                    logService.complete(taskId);

                } catch (Exception e) {
                    logService.sendLog(taskId, "错误: " + e.getMessage());
                }

            }).start();

            result.put("code", 200);
            result.put("msg", "已启动运行");

        } catch (Exception e) {
            result.put("code", 500);
            result.put("msg", e.getMessage());
        }

        return result;
    }
}