package com.phenglei.util;

import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class HyparaFileEditor {

    private HyparaFileEditor() {
    }

    public static Object parseValue(String raw) {
        String value = raw.trim().replaceAll(";$", "").trim();

        if (value.startsWith("\"") && value.endsWith("\"")) {
            return value.substring(1, value.length() - 1);
        }

        if (value.startsWith("[") && value.endsWith("]")) {
            return value;
        }

        try {
            if (value.contains(".") || value.toLowerCase().contains("e")) {
                return Double.parseDouble(value);
            }
            return Integer.parseInt(value);
        } catch (Exception e) {
            return value;
        }
    }

    public static String formatValue(String type, Object value) {
        return switch (type) {
            case "string" -> {
                String s = String.valueOf(value);
                if (s.startsWith("\"") && s.endsWith("\"")) {
                    yield s;
                }
                yield "\"" + s + "\"";
            }
            case "int" -> String.valueOf(toInteger(value));
            case "double" -> String.valueOf(toDouble(value));
            case "array" -> formatArrayValue(value);
            default -> throw new IllegalArgumentException("Unsupported type: " + type);
        };
    }

    private static Integer toInteger(Object value) {
        if (value instanceof Integer i) return i;
        if (value instanceof Long l) return l.intValue();
        if (value instanceof Double d) return d.intValue();
        if (value instanceof Float f) return f.intValue();
        if (value instanceof String s) return Integer.parseInt(s.trim());
        throw new IllegalArgumentException("Cannot convert to int: " + value);
    }

    private static Double toDouble(Object value) {
        if (value instanceof Double d) return d;
        if (value instanceof Float f) return (double) f;
        if (value instanceof Integer i) return i.doubleValue();
        if (value instanceof Long l) return l.doubleValue();
        if (value instanceof String s) return Double.parseDouble(s.trim());
        throw new IllegalArgumentException("Cannot convert to double: " + value);
    }

    private static String formatArrayValue(Object value) {
        if (value instanceof List<?> list) {
            List<String> items = new ArrayList<>();
            for (Object item : list) {
                items.add(String.valueOf(item));
            }
            return "[" + String.join(", ", items) + "]";
        }
        String s = String.valueOf(value).trim();
        if (s.startsWith("[") && s.endsWith("]")) {
            return s;
        }
        return "[" + s + "]";
    }

    public static Map<String, Map<String, Object>> parseAllAssignmentsWithType(String text) {
        Map<String, Map<String, Object>> result = new LinkedHashMap<>();

        Pattern pattern = Pattern.compile(
                "(?m)^\\s*(int|double|string)\\s+([A-Za-z_][A-Za-z0-9_]*)\\s*(\\[\\])?\\s*=\\s*(.*?)\\s*;\\s*(?://.*)?$"
        );

        Matcher matcher = pattern.matcher(text);
        while (matcher.find()) {
            String declaredType = matcher.group(1);
            String key = matcher.group(2);
            String isArray = matcher.group(3);
            String rawValue = matcher.group(4);

            Map<String, Object> item = new LinkedHashMap<>();
            if (isArray != null) {
                item.put("type", "array");
                item.put("value", rawValue.trim());
            } else {
                item.put("type", declaredType);
                item.put("value", parseValue(rawValue));
            }
            result.put(key, item);
        }

        return result;
    }

    public static Map<String, Object> filterValues(Map<String, Map<String, Object>> source, Set<String> keys) {
        Map<String, Object> result = new LinkedHashMap<>();
        for (String key : keys) {
            if (source.containsKey(key)) {
                result.put(key, source.get(key).get("value"));
            }
        }
        return result;
    }

    public static Map<String, Object> parseBoundaryFile(
            String text,
            Set<String> boundaryGlobalKeys,
            Set<String> boundaryBlockKeys
    ) {
        Map<String, Object> result = new LinkedHashMap<>();

        Map<String, Map<String, Object>> globalParsed = parseAllAssignmentsWithType(text);
        for (String key : boundaryGlobalKeys) {
            if (globalParsed.containsKey(key)) {
                result.put(key, globalParsed.get(key).get("value"));
            }
        }

        Map<String, Object> boundaries = new LinkedHashMap<>();

        Pattern blockPattern = Pattern.compile(
                "string\\s+bcName\\s*=\\s*\"([^\"]+)\"\\s*;\\s*\\{(.*?)\\n\\}",
                Pattern.DOTALL
        );
        Matcher blockMatcher = blockPattern.matcher(text);

        while (blockMatcher.find()) {
            String bcName = blockMatcher.group(1);
            String body = blockMatcher.group(2);

            Map<String, Map<String, Object>> parsedBlock = parseAllAssignmentsWithType(body);
            Map<String, Object> filtered = new LinkedHashMap<>();
            for (String key : boundaryBlockKeys) {
                if (parsedBlock.containsKey(key)) {
                    filtered.put(key, parsedBlock.get(key).get("value"));
                }
            }
            boundaries.put(bcName, filtered);
        }

        result.put("boundaries", boundaries);
        return result;
    }

    public static String replaceGlobalAssignment(String text, String key, Object value) {
        String currentType = findType(text, key);
        if (currentType == null) {
            throw new IllegalArgumentException("Parameter not found: " + key);
        }

        String newValue = formatValue(currentType, value);

        Pattern pattern;
        if ("array".equals(currentType)) {
            pattern = Pattern.compile(
                    "(?m)^(\\s*(?:int|double)\\s+" + Pattern.quote(key) + "\\[\\]\\s*=\\s*)(.*?)(\\s*;\\s*(?://.*)?$)"
            );
        } else {
            pattern = Pattern.compile(
                    "(?m)^(\\s*(?:int|double|string)\\s+" + Pattern.quote(key) + "\\s*=\\s*)(.*?)(\\s*;\\s*(?://.*)?$)"
            );
        }

        Matcher matcher = pattern.matcher(text);
        if (!matcher.find()) {
            throw new IllegalArgumentException("Parameter not found: " + key);
        }

        return matcher.replaceFirst(Matcher.quoteReplacement(
                matcher.group(1) + newValue + matcher.group(3)
        ));
    }

    public static String replaceBoundaryAssignment(String text, String bcName, String key, Object value) {
        Pattern blockPattern = Pattern.compile(
                "string\\s+bcName\\s*=\\s*\"" + Pattern.quote(bcName) + "\"\\s*;\\s*\\{(.*?)\\n\\}",
                Pattern.DOTALL
        );
        Matcher blockMatcher = blockPattern.matcher(text);

        if (!blockMatcher.find()) {
            throw new IllegalArgumentException("Boundary block not found: " + bcName);
        }

        String body = blockMatcher.group(1);
        String currentType = findType(body, key);
        if (currentType == null) {
            throw new IllegalArgumentException("Parameter not found in boundary block: " + key);
        }

        String newValue = formatValue(currentType, value);

        Pattern valuePattern;
        if ("array".equals(currentType)) {
            valuePattern = Pattern.compile(
                    "(?m)^(\\s*(?:int|double)\\s+" + Pattern.quote(key) + "\\[\\]\\s*=\\s*)(.*?)(\\s*;\\s*(?://.*)?$)"
            );
        } else {
            valuePattern = Pattern.compile(
                    "(?m)^(\\s*(?:int|double|string)\\s+" + Pattern.quote(key) + "\\s*=\\s*)(.*?)(\\s*;\\s*(?://.*)?$)"
            );
        }

        Matcher valueMatcher = valuePattern.matcher(body);
        if (!valueMatcher.find()) {
            throw new IllegalArgumentException("Parameter not found in boundary block: " + key);
        }

        String newBody = valueMatcher.replaceFirst(Matcher.quoteReplacement(
                valueMatcher.group(1) + newValue + valueMatcher.group(3)
        ));

        return text.substring(0, blockMatcher.start(1)) + newBody + text.substring(blockMatcher.end(1));
    }

    public static String findType(String text, String key) {
        Pattern arrayPattern = Pattern.compile(
                "(?m)^\\s*(int|double)\\s+" + Pattern.quote(key) + "\\[\\]\\s*=\\s*.*?\\s*;\\s*(?://.*)?$"
        );
        Matcher arrayMatcher = arrayPattern.matcher(text);
        if (arrayMatcher.find()) {
            return "array";
        }

        Pattern scalarPattern = Pattern.compile(
                "(?m)^\\s*(int|double|string)\\s+" + Pattern.quote(key) + "\\s*=\\s*.*?\\s*;\\s*(?://.*)?$"
        );
        Matcher scalarMatcher = scalarPattern.matcher(text);
        if (scalarMatcher.find()) {
            return scalarMatcher.group(1);
        }

        return null;
    }
}