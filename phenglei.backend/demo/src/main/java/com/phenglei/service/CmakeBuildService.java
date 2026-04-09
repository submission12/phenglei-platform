package com.phenglei.service;

import com.phenglei.util.FileUtil;
import lombok.RequiredArgsConstructor;
import lombok.extern.slf4j.Slf4j;
import org.springframework.beans.factory.annotation.Value;
import org.springframework.scheduling.annotation.Async;
import org.springframework.stereotype.Service;

import java.io.BufferedReader;
import java.io.File;
import java.io.InputStreamReader;
import java.nio.charset.Charset;
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

    /**
     * 编译 + 拷贝 EXE + DLL
     */
    @Async("taskExecutor")
    public CompletableFuture<Boolean> buildExe() {
        try {
            // 初始化目录
            FileUtil.deleteDir(buildDir);
            FileUtil.createDir(buildDir);
            FileUtil.createDir(installDir);
            log.info("初始化完成 buildDir={}, installDir={}", buildDir, installDir);

            // 生成工程
            String generateCmd = String.format(
                    "cmake -S \"%s\" -B \"%s\" -A x64",
                    sourceDir, buildDir
            );
            if (!executeCmd(generateCmd)) return fail("CMake生成失败");

            // 编译
            String buildCmd = String.format(
                    "cmake --build \"%s\" --config Release -- /m:4",
                    buildDir
            );
            if (!executeCmd(buildCmd)) return fail("CMake编译失败");

            // 复制 EXE
            String srcExe = buildDir + "\\PHengLEIv3d0\\Release\\PHengLEIv3d0.exe";
            if (!copyFile(srcExe, installDir)) return fail("EXE复制失败");

            // 复制 DLL
            String srcDll = buildDir + "\\PHengLEIv3d0\\Release\\tecio.dll";
            tecioDllPath = installDir + "\\tecio.dll";
            if (!copyFile(srcDll, installDir)) return fail("DLL复制失败");

            // 校验
            if (!FileUtil.checkFileExists(exePath)) return fail("EXE不存在");
            if (!FileUtil.checkFileExists(tecioDllPath)) return fail("DLL不存在");

            log.info("构建完成 EXE={}, DLL={}", exePath, tecioDllPath);
            return CompletableFuture.completedFuture(true);

        } catch (Exception e) {
            log.error("构建异常", e);
            return CompletableFuture.completedFuture(false);
        }
    }

    /**
     * ✅ 运行 EXE（关键新增）
     */
    @Async("taskExecutor")
    public CompletableFuture<Boolean> runExe() {
        try {
            File exeFile = new File(exePath);

            if (!exeFile.exists()) {
                log.error("EXE不存在: {}", exePath);
                return CompletableFuture.completedFuture(false);
            }

            log.info("启动EXE: {}", exePath);

            ProcessBuilder pb = new ProcessBuilder(exePath);

            // ⭐必须设置工作目录
            pb.directory(exeFile.getParentFile());

            pb.redirectErrorStream(true);

            Process process = pb.start();

            // 打印输出
            try (BufferedReader reader = new BufferedReader(
                    new InputStreamReader(process.getInputStream(), Charset.forName("GBK")))) {

                String line;
                while ((line = reader.readLine()) != null) {
                    log.info("[EXE] {}", line);
                }
            }

            int code = process.waitFor();
            log.info("EXE结束 exitCode={}", code);

            return CompletableFuture.completedFuture(code == 0);

        } catch (Exception e) {
            log.error("运行失败", e);
            return CompletableFuture.completedFuture(false);
        }
    }

    /**
     * 执行命令
     */
    private boolean executeCmd(String cmd) {
        try {
            log.info("执行: {}", cmd);

            Process process = new ProcessBuilder("cmd.exe", "/c", cmd)
                    .redirectErrorStream(true)
                    .start();

            try (BufferedReader reader = new BufferedReader(
                    new InputStreamReader(process.getInputStream(), Charset.forName("GBK")))) {
                String line;
                while ((line = reader.readLine()) != null) {
                    log.info(line);
                }
            }

            return process.waitFor() == 0;

        } catch (Exception e) {
            log.error("命令执行失败", e);
            return false;
        }
    }

    /**
     * 复制文件（CMD + PowerShell 双保险）
     */
    private boolean copyFile(String src, String destDir) {
        String cmd = String.format("cmd.exe /c copy /Y \"%s\" \"%s\"", src, destDir);
        if (executeCmd(cmd)) return true;

        String ps = String.format(
                "powershell.exe -Command Copy-Item -Path \"%s\" -Destination \"%s\" -Force",
                src, destDir
        );
        return executeCmd(ps);
    }

    private CompletableFuture<Boolean> fail(String msg) {
        log.error(msg);
        return CompletableFuture.completedFuture(false);
    }

    public String getExePath() {
        return exePath;
    }

    public String getTecioDllPath() {
        return tecioDllPath;
    }
}