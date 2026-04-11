<template>
  <div class="container">
    <h2>配置管理</h2>

    <!-- ================= KEY ================= -->
    <div class="section">
      <h2>Key 配置</h2>

      <div class="key-grid">

        <div class="key-line">
          <label>ndim:</label>
          <input v-model="key.ndim" :placeholder="keyOrigin.ndim" />
        </div>

        <div class="key-line">
          <label>nparafile:</label>
          <input v-model="key.nparafile" :placeholder="keyOrigin.nparafile" />
        </div>

        <div class="key-line">
          <label>nsimutask:</label>
          <input v-model="key.nsimutask" :placeholder="keyOrigin.nsimutask" />
        </div>

        <div class="key-line">
          <label>parafilename:</label>
          <select v-model="key.parafilename">
            <option disabled value="">
              {{ keyOrigin.parafilename || '请选择' }}
            </option>
            <option v-for="item in paraOptions" :key="item" :value="item">
              {{ item }}
            </option>
          </select>
        </div>

      </div>
    </div>

    <!-- ================= CFD ================= -->
    <div class="section">
      <h2>CFD 配置</h2>

      <div class="cfd-grid">

        <div
            class="param-line"
            v-for="(value, keyName) in cfdOrigin"
            :key="keyName"
        >
          <label>{{ keyName }}:</label>

          <input
              v-model="cfd[keyName]"
              :placeholder="value"
          />
        </div>

      </div>
    </div>

    <!-- ================= BOUNDARY ================= -->
    <div class="section">
      <h2>Boundary 配置</h2>

      <div class="boundary-container">
        <textarea v-model="boundaryText" class="boundary-textarea"></textarea>
        <textarea v-model="boundaryFileText" class="boundary-textarea readonly" readonly></textarea>
      </div>
    </div>

    <!-- ================= BUTTON ================= -->
    <div class="button-bar">
      <button @click="saveKey">保存 Key</button>
      <button @click="saveCfd">保存 CFD</button>
      <button @click="saveBoundary">保存 Boundary</button>
      <button @click="saveAll">保存全部</button>
      <button @click="resetInitialConfig">恢复默认配置</button>
    </div>

  </div>
</template>

<script>
import { ref, onMounted } from 'vue'
import axios from 'axios'
import { saveAllConfig } from '../api/config.js'

export default {
  setup() {

    const key = ref({})
    const keyOrigin = ref({})

    const cfd = ref({})
    const cfdOrigin = ref({})

    const boundaryText = ref('')
    const boundaryFileText = ref('')

    const paraOptions = [
      "./bin/cfd_para_subsonic.hypara",
      "./bin/cfd_para_transonic.hypara",
      "./bin/cfd_para_supersonic.hypara",
      "./bin/cfd_para_hypersonic.hypara",
      "./bin/cfd_para_incompressible.hypara"
    ]

    const loadConfig = async () => {
      const keyRes = await axios.get('http://localhost:8080/config/key')
      const cfdRes = await axios.get('http://localhost:8080/config/cfd')
      const boundaryRes = await axios.get('http://localhost:8080/config/boundary')

      keyOrigin.value = keyRes.data || {}

      key.value = {
        ndim: '',
        nparafile: '',
        nsimutask: '',
        parafilename: keyRes.data?.parafilename || ''
      }

      let cfdData = cfdRes.data
      if (typeof cfdData === 'string') {
        try { cfdData = JSON.parse(cfdData) } catch { cfdData = {} }
      }

      cfdOrigin.value = cfdData || {}

      const empty = {}
      for (const k in cfdOrigin.value) {
        empty[k] = ''
      }
      cfd.value = empty

      const b = typeof boundaryRes.data === 'string'
          ? boundaryRes.data
          : JSON.stringify(boundaryRes.data, null, 2)

      boundaryText.value = b
      boundaryFileText.value = b
    }

    onMounted(loadConfig)

    const mergeConfig = (input, origin) => {
      const result = {}
      for (const k in origin) {
        result[k] = input[k] !== '' ? input[k] : origin[k]
      }
      return result
    }

    const saveKey = async () => {
      await axios.patch('http://localhost:8080/config/key', key.value)
      alert('Key保存成功')
      window.location.reload()
    }

    const saveCfd = async () => {
      const payload = mergeConfig(cfd.value, cfdOrigin.value)
      await axios.patch('http://localhost:8080/config/cfd', payload)
      alert('CFD保存成功')
      window.location.reload()
    }

    const saveBoundary = async () => {
      const obj = JSON.parse(boundaryText.value)
      await axios.patch('http://localhost:8080/config/boundary', obj)
      alert('Boundary保存成功')
      window.location.reload()
    }

    const saveAll = async () => {
      const obj = JSON.parse(boundaryText.value)
      await axios.patch('http://localhost:8080/config/boundary', obj)
      const payload = mergeConfig(cfd.value, cfdOrigin.value)
      await axios.patch('http://localhost:8080/config/cfd', payload)
      await axios.patch('http://localhost:8080/config/key', key.value)

      alert('全部保存成功')
      window.location.reload()
    }
    const resetInitialConfig = async () => {
      const res = await axios.post('http://localhost:8080/config/reset')
      if (res.data?.code === true) {
        alert('恢复成功')
        await loadConfig()
      } else {
        alert('恢复失败')

      }
      window.location.reload()
    }

    return {
      key,
      keyOrigin,
      cfd,
      cfdOrigin,
      boundaryText,
      boundaryFileText,
      paraOptions,

      saveKey,
      saveCfd,
      saveBoundary,
      saveAll,
      resetInitialConfig
    }
  }
}
</script>

<style>
.container {
  padding: 16px;
  font-family: Arial;
}

.section {
  border: 1px solid #ccc;
  padding: 12px;
  margin-bottom: 16px;
}

/* ================= KEY（已改2列） ================= */
.key-grid {
  display: grid;
  grid-template-columns: repeat(2, 1fr);  /* ⭐ 2列 */
  gap: 10px 40px;
}

.key-line {
  display: grid;
  grid-template-columns: 220px 1fr; /* ⭐ 和CFD统一 */
  align-items: center;
  gap: 10px;
}

.key-line label {
  text-align: right;
  font-weight: 500;
  white-space: nowrap;
}

.key-line input,
.key-line select {
  width: 100%;
  height: 28px;
  font-size: 13px;
  padding: 2px 6px;
  box-sizing: border-box;
}

/* ================= CFD（不动） ================= */
.cfd-grid {
  display: grid;
  grid-template-columns: repeat(2, 1fr);
  gap: 10px 40px;
}

.param-line {
  display: grid;
  grid-template-columns: 220px 1fr;
  align-items: center;
  gap: 10px;
}

.param-line label {
  text-align: right;
  font-weight: 500;
  white-space: nowrap;
}

.param-line input,
.param-line select {
  width: 100%;
  box-sizing: border-box;
}

/* ================= Boundary ================= */
.boundary-container {
  display: flex;
  gap: 12px;
}

.boundary-textarea {
  width: 50%;
  height: 350px;
  font-family: monospace;
  padding: 10px;
}

.readonly {
  background: #f5f5f5;
}

/* ================= button ================= */
.button-bar {
  display: flex;
  justify-content: flex-end;
  gap: 20px;
  margin-top: 20px;
}

</style>