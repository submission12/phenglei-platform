import { createRouter, createWebHistory } from 'vue-router'
import ConfigView from '../views/ConfigView.vue'
import BuildView from '../views/BuildView.vue'

const routes = [
    { path: '/', name: 'Config', component: ConfigView },
    { path: '/build', name: 'Build', component: BuildView }
]

const router = createRouter({
    history: createWebHistory(),
    routes,
})

export default router