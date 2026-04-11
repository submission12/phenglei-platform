<template>
  <div>
    <textarea v-model="jsonText" rows="20" style="width:100%;"></textarea>
  </div>
</template>

<script>
export default {
  name: 'JsonEditor',
  props: {
    value: {
      type: Object,
      required: true
    }
  },
  data() {
    return {
      jsonText: JSON.stringify(this.value, null, 2)
    }
  },
  watch: {
    value: {
      deep: true,
      handler(newVal) {
        this.jsonText = JSON.stringify(newVal, null, 2)
      }
    },
    jsonText(newVal) {
      try {
        this.$emit('update:value', JSON.parse(newVal))
      } catch (e) {
        // 无效 JSON，不触发更新
      }
    }
  }
}
</script>