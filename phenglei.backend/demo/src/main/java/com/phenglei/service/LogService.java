package com.phenglei.service;

import lombok.extern.slf4j.Slf4j;
import org.springframework.stereotype.Service;
import org.springframework.web.servlet.mvc.method.annotation.SseEmitter;

import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;

@Slf4j
@Service
public class LogService {

    private final Map<String, SseEmitter> emitters = new ConcurrentHashMap<>();

    public SseEmitter connect(String taskId) {

        SseEmitter emitter = new SseEmitter(0L);

        emitters.put(taskId, emitter);

        emitter.onCompletion(() -> emitters.remove(taskId));
        emitter.onTimeout(() -> emitters.remove(taskId));
        emitter.onError(e -> emitters.remove(taskId));

        return emitter;
    }

    public void sendLog(String taskId, String msg) {
        try {
            SseEmitter emitter = emitters.get(taskId);
            if (emitter != null) {
                emitter.send(SseEmitter.event().name("log").data(msg));
            }
        } catch (Exception e) {
            emitters.remove(taskId);
        }
    }

    public void sendProgress(String taskId, int progress) {
        try {
            SseEmitter emitter = emitters.get(taskId);
            if (emitter != null) {
                emitter.send(SseEmitter.event().name("progress").data(progress));
            }
        } catch (Exception e) {
            emitters.remove(taskId);
        }
    }

    public void complete(String taskId) {
        SseEmitter emitter = emitters.remove(taskId);
        if (emitter != null) {
            emitter.complete();
        }
    }
}