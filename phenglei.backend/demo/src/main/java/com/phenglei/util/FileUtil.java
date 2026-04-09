package com.phenglei.util;

import lombok.extern.slf4j.Slf4j;
import org.springframework.stereotype.Component;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;

@Slf4j
@Component
public class FileUtil {

    /**
     * 删除目录（包含所有子文件/子目录）
     */
    public static boolean deleteDir(String dirPath) {
        File dir = new File(dirPath);
        if (!dir.exists()) {
            return true;
        }
        File[] files = dir.listFiles();
        if (files != null) {
            for (File file : files) {
                if (file.isDirectory()) {
                    deleteDir(file.getAbsolutePath());
                } else {
                    if (!file.delete()) {
                        log.warn("删除文件失败：{}", file.getAbsolutePath());
                    }
                }
            }
        }
        return dir.delete();
    }

    /**
     * 创建目录（递归创建父目录）
     */
    public static boolean createDir(String dirPath) {
        Path path = Paths.get(dirPath);
        try {
            Files.createDirectories(path);
            return true;
        } catch (Exception e) {
            log.error("创建目录失败：{}", dirPath, e);
            return false;
        }
    }

    /**
     * 检查文件是否存在（且是文件，不是目录）
     */
    public static boolean checkFileExists(String filePath) {
        File file = new File(filePath);
        return file.exists() && file.isFile();
    }
}