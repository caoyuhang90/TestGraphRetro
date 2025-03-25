import { createRouter, createWebHashHistory } from "vue-router";

//路由配置  路由就是用来管理网页之间跳转的机制
const routes = [
    {
        path: '/',
        redirect: '/HomePage'//重定向
    },

    {
        // URL 路径   http://xxxx.com/login 时，这个路由会被匹配
        path: '/login',
        // 定义名字，方便在代码中引用
        name: 'login',


        // 动态导入可以在运行时加载模块  ；  在用户访问该路由时，才去加载 index.vue 组件
        // component 属性指定了该路由对应的组件
        component: () => import('../views/login/index.vue')
    },


    {
        path: '/home',
        name: 'home',
        component: () => import('../views/home/Index.vue'),



        // 嵌套路由
        children: [
            {
                path: '/test0',
                name: 'test0',
                component: () => import('../views/pro0/index.vue'),
                meta: { title: 'test0' }
            },

            {
                path: '/ZiXuanDuanJian',
                name: 'ZiXuanDuanJian',
                component: () => import('../views/ZiXuanDuanJian/index.vue'),
                meta: { title: '自选断键' }
            },

            {
                path: '/HomePage',
                name: 'HomePage',
                component: () => import('../views/home/index/index.vue'),
                meta: { title: '首页' }
            },

            {
                path: '/ChanlvYuce',
                name: 'ChanlvYuce',
                component: () => import('../views/ChanlvYuce/index.vue'),
                meta: { title: ' 产率预测' }
            },

            {
                path: '/ChanwuYuce',
                name: 'ChanwuYuce',
                component: () => import('../views/ChanwuYuce/index.vue'),
                meta: { title: '产物预测' }
            },

            {
                path: '/LeixingYuce',
                name: 'LeixingYuce',
                component: () => import('../views/LeixingYuce/index.vue'),
                meta: { title: '类型预测' }
            },

            {
                path: '/LujingYuce',
                name: 'LujingYuce',
                component: () => import('../views/LujingYuce/index.vue'),
                meta: { title: '路径预测' }
            },

            {
                path: '/TiaojianYouhua',
                name: 'TiaojianYouhua',
                component: () => import('../views/TiaojianYouhua/index.vue'),
                meta: { title: '条件预测' }
            },

            {
                path: '/XinxiChouqu',
                name: 'XinxiChouqu',
                component: () => import('../views/XinxiChouqu/index.vue'),
                meta: { title: '信息抽取' }
            },

        ]
    },

]


const router = createRouter({
    history: createWebHashHistory(),
    routes
})

// 导出 router 实例, 供应用程序其他部分使用
export default router