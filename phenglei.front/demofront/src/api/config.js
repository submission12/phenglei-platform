import axios from 'axios'

const api = axios.create({
    baseURL: '/', // 修改为你后端实际地址
    timeout: 10000
})

export const getKeyConfig = () => api.get('/config/key')
export const getCfdConfig = () => api.get('/config/cfd')
export const getBoundaryConfig = () => api.get('/config/boundary')
export const saveAllConfig = (data) => api.post('/config/save', data)