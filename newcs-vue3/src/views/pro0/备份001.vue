<template>
    <el-container style="height: 100vh;">
        <!-- 顶部导航栏 -->
        <el-header class="header">
            测试0
        </el-header>
        <div>
            <h1>分子结构生成器</h1>
            <input v-model="moleculeInput" placeholder="输入分子表达式" />
            <button @click="jumpTest">跳转测试</button>
            <hr>

            <hr>
            O=C1C=CC2=C(CC[C@H](CC)[C@@H]2O)O1<br>
            O=C1C2([H])C=C(CO)[C@@H](O)[C@]3(O)[C@@H](O)C(C)=C[C@@]31[C@H](C)C[C@]4([H])[C@@]2([H])C4(C)C<br>
            <hr>
            全体图片<br>
            <img v-for="(imgURL, index) in list" :key="index" :src="imgURL" alt="分子结构图" />

            <hr>


        </div>


    </el-container>
</template>


<script setup >
    import { ref } from 'vue';
    import axios from 'axios';

    const moleculeInput = ref('');
    const imageUrl = ref('');
    let result = ref()
    let list = ref([])
    let i = 0
    let imgBase64 = 0

    async function jumpTest() {  // async异步要求
        console.log('1111')

        result = await axios.post('http://127.0.0.1:3000/api/generate-molecule', { molecule: String(moleculeInput.value) })
        if (result.data.images && result.data.images.length > 0) {
            list.value = [];
            for (i = 0; i < 5 || i < result.data.images.length; i++) {
                imgBase64 = `data:image/png;base64,${result.data.images[i]}`;
                list.value.push(imgBase64)
            }
        }

        // 可以执行
        // console.log('响应数据:', result.data)
        // imageUrl.value = `data:image/png;base64,${result.data.images[0]}`;
    }


</script>


<style scoped>
    .header {
        background-color: rgb(1, 70, 138);
        color: white;
        text-align: center;
        padding: 0px 0;
    }

    .el-tabs {
        background-color: #ffffff;
        margin-bottom: 20px;
    }

    .module {
        padding: 20px;
        box-shadow: 0 2px 8px rgba(0, 0, 0, 0.1);
        border-radius: 8px;
    }

    .el-upload__text {
        font-size: 12px;
        text-align: center;
    }

    .el-button {
        width: 100%;
    }
</style>