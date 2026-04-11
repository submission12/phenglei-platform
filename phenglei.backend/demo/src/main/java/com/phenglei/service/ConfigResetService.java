package com.phenglei.service;

import lombok.extern.slf4j.Slf4j;
import org.springframework.beans.factory.annotation.Value;
import org.springframework.stereotype.Service;

import java.io.File;
import java.io.InputStream;
import java.io.OutputStream;
import java.io.FileInputStream;
import java.io.FileOutputStream;

@Service
@Slf4j
public class ConfigResetService {

    @Value("${phenglei.config.default-dir}")
    private String defaultDir;

    @Value("${phenglei.config.key-file}")
    private String keyFile;

    @Value("${phenglei.config.cfd-file}")
    private String cfdFile;

    @Value("${phenglei.config.boundary-file}")
    private String boundaryFile;

    public void resetToDefault() {
        try {
            // 默认文件
            File defaultKey = new File(defaultDir, "key.hypara");
            File defaultCfd = new File(defaultDir, "cfd_para.hypara");
            File defaultBoundary = new File(defaultDir, "boundary_condition.hypara");

            // 目标文件
            File targetKey = new File(keyFile);
            File targetCfd = new File(cfdFile);
            File targetBoundary = new File(boundaryFile);

            // 覆盖复制
            copyFile(defaultKey, targetKey);
            copyFile(defaultCfd, targetCfd);
            copyFile(defaultBoundary, targetBoundary);

            log.info("配置已恢复为默认");

        } catch (Exception e) {
            log.error("恢复默认配置失败", e);
            throw new RuntimeException("恢复失败");
        }
    }

    private void copyFile(File source, File target) throws Exception {
        try (InputStream in = new FileInputStream(source);
             OutputStream out = new FileOutputStream(target)) {

            byte[] buffer = new byte[4096];
            int len;
            while ((len = in.read(buffer)) > 0) {
                out.write(buffer, 0, len);
            }
        }
    }
}