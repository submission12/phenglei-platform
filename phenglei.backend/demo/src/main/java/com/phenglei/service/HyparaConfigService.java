package com.phenglei.service;

import com.phenglei.model.BoundaryUpdateRequest;
import com.phenglei.model.CfdUpdateRequest;
import com.phenglei.model.KeyUpdateRequest;
import com.phenglei.util.HyparaFileEditor;
import org.springframework.beans.factory.annotation.Value;
import org.springframework.stereotype.Service;

import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.Map;
import java.util.Set;

@Service
public class HyparaConfigService {

    @Value("${phenglei.config.key-file}")
    private String keyFile;

    @Value("${phenglei.config.cfd-file}")
    private String cfdFile;

    @Value("${phenglei.config.boundary-file}")
    private String boundaryFile;

    private static final Set<String> KEY_KEYS = new LinkedHashSet<>();
    private static final Set<String> CFD_KEYS = new LinkedHashSet<>();
    private static final Set<String> BOUNDARY_GLOBAL_KEYS = new LinkedHashSet<>();
    private static final Set<String> BOUNDARY_BLOCK_KEYS = new LinkedHashSet<>();

    static {
        KEY_KEYS.add("ndim");
        KEY_KEYS.add("nparafile");
        KEY_KEYS.add("nsimutask");
        KEY_KEYS.add("parafilename");

        CFD_KEYS.add("maxSimuStep");
        CFD_KEYS.add("intervalStepFlow");
        CFD_KEYS.add("intervalStepPlot");
        CFD_KEYS.add("intervalStepForce");
        CFD_KEYS.add("intervalStepRes");
        CFD_KEYS.add("ifLowSpeedPrecon");
        CFD_KEYS.add("refMachNumber");
        CFD_KEYS.add("attackd");
        CFD_KEYS.add("angleSlide");
        CFD_KEYS.add("inflowParaType");
        CFD_KEYS.add("refReNumber");
        CFD_KEYS.add("refDimensionalTemperature");
        CFD_KEYS.add("gridScaleFactor");
        CFD_KEYS.add("forceReferenceLengthSpanWise");
        CFD_KEYS.add("forceReferenceLength");
        CFD_KEYS.add("forceReferenceArea");
        CFD_KEYS.add("TorqueRefX");
        CFD_KEYS.add("TorqueRefY");
        CFD_KEYS.add("TorqueRefZ");
        CFD_KEYS.add("viscousType");
        CFD_KEYS.add("viscousName");
        CFD_KEYS.add("DESType");
        CFD_KEYS.add("roeEntropyFixMethod");
        CFD_KEYS.add("roeEntropyScale");
        CFD_KEYS.add("str_limiter_name");
        CFD_KEYS.add("uns_limiter_name");
        CFD_KEYS.add("venkatCoeff");
        CFD_KEYS.add("iunsteady");
        CFD_KEYS.add("CFLEnd");
        CFD_KEYS.add("nLUSGSSweeps");
        CFD_KEYS.add("nMGLevel");
        CFD_KEYS.add("flowInitStep");
        CFD_KEYS.add("gridfile");
        CFD_KEYS.add("plotFieldType");
        CFD_KEYS.add("nVisualVariables");
        CFD_KEYS.add("visualVariables");
        CFD_KEYS.add("reconmeth");
        CFD_KEYS.add("limitVariables");
        CFD_KEYS.add("limitVector");

        BOUNDARY_GLOBAL_KEYS.add("nBoundaryConditions");

        BOUNDARY_BLOCK_KEYS.add("bcType");
        BOUNDARY_BLOCK_KEYS.add("forceReferenceLength");
        BOUNDARY_BLOCK_KEYS.add("forceReferenceLengthSpanWise");
        BOUNDARY_BLOCK_KEYS.add("forceReferenceArea");
        BOUNDARY_BLOCK_KEYS.add("TorqueRefX");
        BOUNDARY_BLOCK_KEYS.add("TorqueRefY");
        BOUNDARY_BLOCK_KEYS.add("TorqueRefZ");
        BOUNDARY_BLOCK_KEYS.add("dumpHingeMoment");
        BOUNDARY_BLOCK_KEYS.add("localCoordAxis0");
        BOUNDARY_BLOCK_KEYS.add("localCoordAxis1");

        BOUNDARY_BLOCK_KEYS.add("viscousType");
        BOUNDARY_BLOCK_KEYS.add("wallTemperature");
        BOUNDARY_BLOCK_KEYS.add("uWall");
        BOUNDARY_BLOCK_KEYS.add("vWall");
        BOUNDARY_BLOCK_KEYS.add("wWall");

        BOUNDARY_BLOCK_KEYS.add("inflowParaType");
        BOUNDARY_BLOCK_KEYS.add("refMachNumber");
        BOUNDARY_BLOCK_KEYS.add("attackd");
        BOUNDARY_BLOCK_KEYS.add("angleSlide");
        BOUNDARY_BLOCK_KEYS.add("refReNumber");
        BOUNDARY_BLOCK_KEYS.add("refDimensionalTemperature");
        BOUNDARY_BLOCK_KEYS.add("refDimensionalPressure");
        BOUNDARY_BLOCK_KEYS.add("refDimensionalVelocity");
        BOUNDARY_BLOCK_KEYS.add("refDimensionalDensity");
        BOUNDARY_BLOCK_KEYS.add("refDimensionalHeight");
        BOUNDARY_BLOCK_KEYS.add("powerLawCoefficient");
        BOUNDARY_BLOCK_KEYS.add("height");

        BOUNDARY_BLOCK_KEYS.add("primDensity");
        BOUNDARY_BLOCK_KEYS.add("primU");
        BOUNDARY_BLOCK_KEYS.add("primV");
        BOUNDARY_BLOCK_KEYS.add("primW");
        BOUNDARY_BLOCK_KEYS.add("primPressure");

        BOUNDARY_BLOCK_KEYS.add("nTrajectoryVariables");
        BOUNDARY_BLOCK_KEYS.add("time");
        BOUNDARY_BLOCK_KEYS.add("pressure");
        BOUNDARY_BLOCK_KEYS.add("velocity");
        BOUNDARY_BLOCK_KEYS.add("temperature");

        BOUNDARY_BLOCK_KEYS.add("totalPressure");
        BOUNDARY_BLOCK_KEYS.add("totalTemperature");
        BOUNDARY_BLOCK_KEYS.add("direction_inlet");
        BOUNDARY_BLOCK_KEYS.add("staticPressure");
        BOUNDARY_BLOCK_KEYS.add("massFlow");
    }

    public Map<String, Object> getKeyConfig() throws IOException {
        String text = Files.readString(Path.of(keyFile), StandardCharsets.UTF_8);
        return HyparaFileEditor.filterValues(
                HyparaFileEditor.parseAllAssignmentsWithType(text),
                KEY_KEYS
        );
    }

    public Map<String, Object> getCfdConfig() throws IOException {
        String text = Files.readString(Path.of(cfdFile), StandardCharsets.UTF_8);
        return HyparaFileEditor.filterValues(
                HyparaFileEditor.parseAllAssignmentsWithType(text),
                CFD_KEYS
        );
    }

    public Map<String, Object> getBoundaryConfig() throws IOException {
        String text = Files.readString(Path.of(boundaryFile), StandardCharsets.UTF_8);
        return HyparaFileEditor.parseBoundaryFile(text, BOUNDARY_GLOBAL_KEYS, BOUNDARY_BLOCK_KEYS);
    }

    public Map<String, Object> updateKeyConfig(KeyUpdateRequest request) throws IOException {
        Path path = Path.of(keyFile);
        String text = Files.readString(path, StandardCharsets.UTF_8);

        boolean changed = false;
        Map<String, Object> updatedValues = new LinkedHashMap<>();

        if (request.getNdim() != null) {
            if (HyparaFileEditor.findType(text, "ndim") == null) {
                throw new IllegalArgumentException("key.hypara 中未找到参数: ndim");
            }
            text = HyparaFileEditor.replaceGlobalAssignment(text, "ndim", request.getNdim());
            updatedValues.put("ndim", request.getNdim());
            changed = true;
        }

        if (request.getNparafile() != null) {
            if (HyparaFileEditor.findType(text, "nparafile") == null) {
                throw new IllegalArgumentException("key.hypara 中未找到参数: nparafile");
            }
            text = HyparaFileEditor.replaceGlobalAssignment(text, "nparafile", request.getNparafile());
            updatedValues.put("nparafile", request.getNparafile());
            changed = true;
        }

        if (request.getNsimutask() != null) {
            if (HyparaFileEditor.findType(text, "nsimutask") == null) {
                throw new IllegalArgumentException("key.hypara 中未找到参数: nsimutask");
            }
            text = HyparaFileEditor.replaceGlobalAssignment(text, "nsimutask", request.getNsimutask());
            updatedValues.put("nsimutask", request.getNsimutask());
            changed = true;
        }

        if (request.getParafilename() != null) {
            if (HyparaFileEditor.findType(text, "parafilename") == null) {
                throw new IllegalArgumentException("key.hypara 中未找到参数: parafilename");
            }
            text = HyparaFileEditor.replaceGlobalAssignment(text, "parafilename", request.getParafilename());
            updatedValues.put("parafilename", request.getParafilename());
            changed = true;
        }

        if (!changed) {
            throw new IllegalArgumentException("至少提交一个可修改参数: ndim, nparafile, nsimutask, parafilename");
        }

        Files.writeString(path, text, StandardCharsets.UTF_8);

        Map<String, Object> result = new LinkedHashMap<>();
        result.put("ok", true);
        result.put("file", "key.hypara");
        result.put("updated", updatedValues);
        return result;
    }

    public Map<String, Object> updateCfdConfig(CfdUpdateRequest request) throws IOException {
        Path path = Path.of(cfdFile);
        String text = Files.readString(path, StandardCharsets.UTF_8);

        boolean changed = false;
        Map<String, Object> updatedValues = new LinkedHashMap<>();

        if (request.getMaxSimuStep() != null) {
            text = updateCfdField(text, "maxSimuStep", request.getMaxSimuStep());
            updatedValues.put("maxSimuStep", request.getMaxSimuStep());
            changed = true;
        }
        if (request.getIntervalStepFlow() != null) {
            text = updateCfdField(text, "intervalStepFlow", request.getIntervalStepFlow());
            updatedValues.put("intervalStepFlow", request.getIntervalStepFlow());
            changed = true;
        }
        if (request.getIntervalStepPlot() != null) {
            text = updateCfdField(text, "intervalStepPlot", request.getIntervalStepPlot());
            updatedValues.put("intervalStepPlot", request.getIntervalStepPlot());
            changed = true;
        }
        if (request.getIntervalStepForce() != null) {
            text = updateCfdField(text, "intervalStepForce", request.getIntervalStepForce());
            updatedValues.put("intervalStepForce", request.getIntervalStepForce());
            changed = true;
        }
        if (request.getIntervalStepRes() != null) {
            text = updateCfdField(text, "intervalStepRes", request.getIntervalStepRes());
            updatedValues.put("intervalStepRes", request.getIntervalStepRes());
            changed = true;
        }
        if (request.getIfLowSpeedPrecon() != null) {
            text = updateCfdField(text, "ifLowSpeedPrecon", request.getIfLowSpeedPrecon());
            updatedValues.put("ifLowSpeedPrecon", request.getIfLowSpeedPrecon());
            changed = true;
        }
        if (request.getRefMachNumber() != null) {
            text = updateCfdField(text, "refMachNumber", request.getRefMachNumber());
            updatedValues.put("refMachNumber", request.getRefMachNumber());
            changed = true;
        }
        if (request.getAttackd() != null) {
            text = updateCfdField(text, "attackd", request.getAttackd());
            updatedValues.put("attackd", request.getAttackd());
            changed = true;
        }
        if (request.getAngleSlide() != null) {
            text = updateCfdField(text, "angleSlide", request.getAngleSlide());
            updatedValues.put("angleSlide", request.getAngleSlide());
            changed = true;
        }
        if (request.getInflowParaType() != null) {
            text = updateCfdField(text, "inflowParaType", request.getInflowParaType());
            updatedValues.put("inflowParaType", request.getInflowParaType());
            changed = true;
        }
        if (request.getRefReNumber() != null) {
            text = updateCfdField(text, "refReNumber", request.getRefReNumber());
            updatedValues.put("refReNumber", request.getRefReNumber());
            changed = true;
        }
        if (request.getRefDimensionalTemperature() != null) {
            text = updateCfdField(text, "refDimensionalTemperature", request.getRefDimensionalTemperature());
            updatedValues.put("refDimensionalTemperature", request.getRefDimensionalTemperature());
            changed = true;
        }
        if (request.getGridScaleFactor() != null) {
            text = updateCfdField(text, "gridScaleFactor", request.getGridScaleFactor());
            updatedValues.put("gridScaleFactor", request.getGridScaleFactor());
            changed = true;
        }
        if (request.getForceReferenceLengthSpanWise() != null) {
            text = updateCfdField(text, "forceReferenceLengthSpanWise", request.getForceReferenceLengthSpanWise());
            updatedValues.put("forceReferenceLengthSpanWise", request.getForceReferenceLengthSpanWise());
            changed = true;
        }
        if (request.getForceReferenceLength() != null) {
            text = updateCfdField(text, "forceReferenceLength", request.getForceReferenceLength());
            updatedValues.put("forceReferenceLength", request.getForceReferenceLength());
            changed = true;
        }
        if (request.getForceReferenceArea() != null) {
            text = updateCfdField(text, "forceReferenceArea", request.getForceReferenceArea());
            updatedValues.put("forceReferenceArea", request.getForceReferenceArea());
            changed = true;
        }
        if (request.getTorqueRefX() != null) {
            text = updateCfdField(text, "TorqueRefX", request.getTorqueRefX());
            updatedValues.put("TorqueRefX", request.getTorqueRefX());
            changed = true;
        }
        if (request.getTorqueRefY() != null) {
            text = updateCfdField(text, "TorqueRefY", request.getTorqueRefY());
            updatedValues.put("TorqueRefY", request.getTorqueRefY());
            changed = true;
        }
        if (request.getTorqueRefZ() != null) {
            text = updateCfdField(text, "TorqueRefZ", request.getTorqueRefZ());
            updatedValues.put("TorqueRefZ", request.getTorqueRefZ());
            changed = true;
        }
        if (request.getViscousType() != null) {
            text = updateCfdField(text, "viscousType", request.getViscousType());
            updatedValues.put("viscousType", request.getViscousType());
            changed = true;
        }
        if (request.getViscousName() != null) {
            text = updateCfdField(text, "viscousName", request.getViscousName());
            updatedValues.put("viscousName", request.getViscousName());
            changed = true;
        }
        if (request.getDESType() != null) {
            text = updateCfdField(text, "DESType", request.getDESType());
            updatedValues.put("DESType", request.getDESType());
            changed = true;
        }
        if (request.getRoeEntropyFixMethod() != null) {
            text = updateCfdField(text, "roeEntropyFixMethod", request.getRoeEntropyFixMethod());
            updatedValues.put("roeEntropyFixMethod", request.getRoeEntropyFixMethod());
            changed = true;
        }
        if (request.getRoeEntropyScale() != null) {
            text = updateCfdField(text, "roeEntropyScale", request.getRoeEntropyScale());
            updatedValues.put("roeEntropyScale", request.getRoeEntropyScale());
            changed = true;
        }
        if (request.getStr_limiter_name() != null) {
            text = updateCfdField(text, "str_limiter_name", request.getStr_limiter_name());
            updatedValues.put("str_limiter_name", request.getStr_limiter_name());
            changed = true;
        }
        if (request.getUns_limiter_name() != null) {
            text = updateCfdField(text, "uns_limiter_name", request.getUns_limiter_name());
            updatedValues.put("uns_limiter_name", request.getUns_limiter_name());
            changed = true;
        }
        if (request.getVenkatCoeff() != null) {
            text = updateCfdField(text, "venkatCoeff", request.getVenkatCoeff());
            updatedValues.put("venkatCoeff", request.getVenkatCoeff());
            changed = true;
        }
        if (request.getIunsteady() != null) {
            text = updateCfdField(text, "iunsteady", request.getIunsteady());
            updatedValues.put("iunsteady", request.getIunsteady());
            changed = true;
        }
        if (request.getCFLEnd() != null) {
            text = updateCfdField(text, "CFLEnd", request.getCFLEnd());
            updatedValues.put("CFLEnd", request.getCFLEnd());
            changed = true;
        }
        if (request.getNLUSGSSweeps() != null) {
            text = updateCfdField(text, "nLUSGSSweeps", request.getNLUSGSSweeps());
            updatedValues.put("nLUSGSSweeps", request.getNLUSGSSweeps());
            changed = true;
        }
        if (request.getNMGLevel() != null) {
            text = updateCfdField(text, "nMGLevel", request.getNMGLevel());
            updatedValues.put("nMGLevel", request.getNMGLevel());
            changed = true;
        }
        if (request.getFlowInitStep() != null) {
            text = updateCfdField(text, "flowInitStep", request.getFlowInitStep());
            updatedValues.put("flowInitStep", request.getFlowInitStep());
            changed = true;
        }
        if (request.getGridfile() != null) {
            text = updateCfdField(text, "gridfile", request.getGridfile());
            updatedValues.put("gridfile", request.getGridfile());
            changed = true;
        }
        if (request.getPlotFieldType() != null) {
            text = updateCfdField(text, "plotFieldType", request.getPlotFieldType());
            updatedValues.put("plotFieldType", request.getPlotFieldType());
            changed = true;
        }
        if (request.getNVisualVariables() != null) {
            text = updateCfdField(text, "nVisualVariables", request.getNVisualVariables());
            updatedValues.put("nVisualVariables", request.getNVisualVariables());
            changed = true;
        }
        if (request.getVisualVariables() != null) {
            text = updateCfdField(text, "visualVariables", request.getVisualVariables());
            updatedValues.put("visualVariables", request.getVisualVariables());
            changed = true;
        }
        if (request.getReconmeth() != null) {
            text = updateCfdField(text, "reconmeth", request.getReconmeth());
            updatedValues.put("reconmeth", request.getReconmeth());
            changed = true;
        }
        if (request.getLimitVariables() != null) {
            text = updateCfdField(text, "limitVariables", request.getLimitVariables());
            updatedValues.put("limitVariables", request.getLimitVariables());
            changed = true;
        }
        if (request.getLimitVector() != null) {
            text = updateCfdField(text, "limitVector", request.getLimitVector());
            updatedValues.put("limitVector", request.getLimitVector());
            changed = true;
        }

        if (!changed) {
            throw new IllegalArgumentException("至少提交一个可修改 CFD 参数");
        }

        Files.writeString(path, text, StandardCharsets.UTF_8);

        Map<String, Object> result = new LinkedHashMap<>();
        result.put("ok", true);
        result.put("file", "cfd_para.hypara");
        result.put("updated", updatedValues);
        return result;
    }

    public Map<String, Object> updateBoundaryConfig(BoundaryUpdateRequest request) throws IOException {
        Path path = Path.of(boundaryFile);
        String text = Files.readString(path, StandardCharsets.UTF_8);

        boolean changed = false;
        Map<String, Object> updated = new LinkedHashMap<>();

        if (request.getNBoundaryConditions() != null) {
            if (!BOUNDARY_GLOBAL_KEYS.contains("nBoundaryConditions")) {
                throw new IllegalArgumentException("boundary_condition.hypara 中不允许修改全局参数: nBoundaryConditions");
            }
            if (HyparaFileEditor.findType(text, "nBoundaryConditions") == null) {
                throw new IllegalArgumentException("boundary_condition.hypara 中未找到参数: nBoundaryConditions");
            }
            text = HyparaFileEditor.replaceGlobalAssignment(text, "nBoundaryConditions", request.getNBoundaryConditions());
            updated.put("nBoundaryConditions", request.getNBoundaryConditions());
            changed = true;
        }

        if (request.getBoundaries() != null) {
            Map<String, Object> updatedBoundaries = new LinkedHashMap<>();

            for (Map.Entry<String, Map<String, Object>> boundaryEntry : request.getBoundaries().entrySet()) {
                String bcName = boundaryEntry.getKey();
                Map<String, Object> fields = boundaryEntry.getValue();

                if (fields == null || fields.isEmpty()) {
                    continue;
                }

                Map<String, Object> updatedFields = new LinkedHashMap<>();

                for (Map.Entry<String, Object> fieldEntry : fields.entrySet()) {
                    String key = fieldEntry.getKey();
                    Object value = fieldEntry.getValue();

                    if (!BOUNDARY_BLOCK_KEYS.contains(key)) {
                        throw new IllegalArgumentException("boundary_condition.hypara 中不允许修改该参数: " + key);
                    }

                    text = HyparaFileEditor.replaceBoundaryAssignment(text, bcName, key, value);
                    updatedFields.put(key, value);
                    changed = true;
                }

                if (!updatedFields.isEmpty()) {
                    updatedBoundaries.put(bcName, updatedFields);
                }
            }

            if (!updatedBoundaries.isEmpty()) {
                updated.put("boundaries", updatedBoundaries);
            }
        }

        if (!changed) {
            throw new IllegalArgumentException("至少提交一个可修改 boundary 参数");
        }

        Files.writeString(path, text, StandardCharsets.UTF_8);

        Map<String, Object> result = new LinkedHashMap<>();
        result.put("ok", true);
        result.put("file", "boundary_condition.hypara");
        result.put("updated", updated);
        return result;
    }

    private String updateCfdField(String text, String key, Object value) {
        if (!CFD_KEYS.contains(key)) {
            throw new IllegalArgumentException("cfd_para.hypara 中不允许修改该参数: " + key);
        }

        if (HyparaFileEditor.findType(text, key) == null) {
            throw new IllegalArgumentException("cfd_para.hypara 中未找到参数: " + key);
        }

        return HyparaFileEditor.replaceGlobalAssignment(text, key, value);
    }
}