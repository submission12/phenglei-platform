package com.phenglei.controller;

import com.phenglei.service.LogService;
import org.springframework.web.bind.annotation.*;
import org.springframework.web.servlet.mvc.method.annotation.SseEmitter;

@RestController
@RequestMapping("/log")
@CrossOrigin
public class LogController {

    private final LogService logService;

    public LogController(LogService logService) {
        this.logService = logService;
    }

    @GetMapping("/subscribe/{taskId}")
    public SseEmitter subscribe(@PathVariable String taskId) {
        return logService.connect(taskId);
    }
}