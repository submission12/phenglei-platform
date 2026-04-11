package com.phenglei.service;

import com.phenglei.util.FileUtil;
import lombok.RequiredArgsConstructor;
import lombok.extern.slf4j.Slf4j;
import org.springframework.beans.factory.annotation.Value;
import org.springframework.scheduling.annotation.Async;
import org.springframework.stereotype.Service;

import java.io.File;
import java.util.concurrent.CompletableFuture;

@Slf4j
@Service
@RequiredArgsConstructor
public class CmakeBuildService {

    @Value("${phenglei.cmake.source-dir}")
    private String sourceDir;

    @Value("${phenglei.cmake.build-dir}")
    private String buildDir;

    @Value("${phenglei.cmake.install-dir}")
    private String installDir;

    @Value("${phenglei.cmake.exe-path}")
    private String exePath;

    private String tecioDllPath;

    @Async("taskExecutor")
    public CompletableFuture<Boolean> buildExe() {
        try {
            FileUtil.deleteDir(buildDir);
            FileUtil.createDir(buildDir);
            FileUtil.createDir(installDir);

            log.info("开始构建");

            exec("cmake -S \"" + sourceDir + "\" -B \"" + buildDir + "\" -A x64");
            exec("cmake --build \"" + buildDir + "\" --config Release -- /m:4");

            return CompletableFuture.completedFuture(true);

        } catch (Exception e) {
            log.error("构建失败", e);
            return CompletableFuture.completedFuture(false);
        }
    }

    /**
     * ⭐ 只负责启动 EXE（关键）
     */
    public Process startExe() throws Exception {

        File exeFile = new File(exePath);

        if (!exeFile.exists()) {
            throw new RuntimeException("EXE不存在: " + exePath);
        }

        ProcessBuilder pb = new ProcessBuilder(exePath);
        pb.directory(exeFile.getParentFile());
        pb.redirectErrorStream(true);

        return pb.start();
    }

    private void exec(String cmd) {
        try {
            Process p = new ProcessBuilder("cmd.exe", "/c", cmd).start();
            p.waitFor();
        } catch (Exception e) {
            throw new RuntimeException(e);
        }
    }

    public String getExePath() {
        return exePath;
    }

    public String getTecioDllPath() {
        return tecioDllPath;
    }
}