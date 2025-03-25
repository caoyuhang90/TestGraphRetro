<template>
    <div class="login">
        <div class="login-content">
            <div class="login-item login-box">
                <i><el-icon><Document /></el-icon></i>
                <p>化学信息抽取</p>
            </div>
            <div class="login-item login-box">
                <i><el-icon><Help /></el-icon></i>
                <p>反应条件优化</p>
            </div>
            <div class="login-form login-item">
                <p>化学合成优化平台</p>
                <el-form
                    ref="ruleFormRef"
                    style="max-width: 600px"
                    :model="ruleForm"
                    :rules="rules"
                    label-width="auto"
                    class="demo-ruleForm"
                >
                    <el-form-item label="user" prop="name">
                    <el-input v-model="ruleForm.name" />
                    </el-form-item>

                    <el-form-item label="password" prop="password">
                    <el-input type="password" v-model="ruleForm.password" />
                    </el-form-item>

                    <el-form-item>
                    <el-button type="success" @click="submitForm(ruleFormRef)">
                        登录
                    </el-button>
                    </el-form-item>
                    
                </el-form>       
            </div>
            <div class="login-item login-box">
                <i><el-icon><MagicStick /></el-icon></i>
                <p>产物产率预测</p>
            </div>
            <div class="login-item login-box">
                <i><el-icon><Guide /></el-icon></i>
                <p>合成路径规划</p>
            </div>

        </div>
        
    </div>
</template>

<script setup>
import {ref,reactive} from 'vue'
import {
    Document,
    Help,
    User,
    MagicStick,
    Guide
} from '@element-plus/icons-vue'
const ruleFormRef = ref();
const ruleForm = reactive({
    name:'',
    password:''
})
const rules = reactive({
    name: [
    { required: true, message: '请输入用户名', trigger: 'blur' },
    { min: 3, max: 12, message: '长度为3-12', trigger: 'blur' },
  ],
    password: [
    { required: true, message: '请输入密码', trigger: 'blur' },
    { min: 8, max: 18, message: '长度为8-18', trigger: 'blur' },
  ],
})
</script>