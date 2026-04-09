package com.phenglei;

import org.springframework.boot.SpringApplication;
import org.springframework.boot.autoconfigure.SpringBootApplication;
import org.springframework.scheduling.annotation.EnableAsync;

@SpringBootApplication

@EnableAsync // 必须加，开启异步任务
public class PhengleiApplication {
    public static void main(String[] args) {
        SpringApplication.run(PhengleiApplication.class, args);
    }
}