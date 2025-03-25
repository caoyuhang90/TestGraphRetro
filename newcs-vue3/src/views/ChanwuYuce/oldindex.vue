<template>
    <div class="container">
      <div class="header">
        <h1>反应产物结构预测</h1>
      </div>
  
      <div class="nav">
        <router-link to="#">创建任务</router-link>
        <router-link to="#">最近结果</router-link>
      </div>
  
      <div class="content">
        <div class="form-container">
          <div class="tab-container">
            <button :class="{ active: activeTab === 'smiles' }" @click="activeTab = 'smiles'">smiles</button>
            <button :class="{ active: activeTab === 'draw' }" @click="activeTab = 'draw'">绘制分子</button>
            <button :class="{ active: activeTab === 'upload' }" @click="activeTab = 'upload'">上传文件</button>
            <button :class="{ active: activeTab === 'data-center' }" @click="activeTab = 'data-center'">数据中心</button>
          </div>
  
          <form @submit.prevent="submitForm">
            <label for="reactants">输入smiles字符串</label>
            <input type="text" v-model="reactants" id="reactants" required />
  
            <label for="task-name">任务名称</label>
            <input type="text" v-model="taskName" id="taskName" />
  
            <label for="reaction-type">反应类型（选填）</label>
            <input type="text" v-model="reactionType" id="reactionType" />
  
            <input type="submit" value="预测产物结构" />
          </form>
        </div>
      </div>
  
      <div class="result-container" v-if="results">
        <h2>预测结果</h2>
        <table>
          <thead>
            <tr>
              <th>Rank</th>
              <th>Product</th>
              <th>Score</th>
              <th>Template</th>
              <th>Pred Actions</th>
            </tr>
          </thead>
          <tbody>
            <tr v-for="(result, key) in results" :key="key">
              <td>{{ key }}</td>
              <td>{{ result.product }}</td>
              <td>{{ result.score }}</td>
              <td>{{ result.template }}</td>
              <td>{{ Array.isArray(result.pred_actions) ? result.pred_actions.join(', ') : '' }}</td>
            </tr>
          </tbody>
        </table>
      </div>
  
      <div class="result-container" v-if="dataFrame">
        <h2>数据框</h2>
        <table>
          <thead>
            <tr>
              <th v-for="header in dataFrame.columns" :key="header">{{ header }}</th>
            </tr>
          </thead>
          <tbody>
            <tr v-for="(row, rowIndex) in dataFrame.data" :key="rowIndex">
              <td v-for="(cell, cellIndex) in row" :key="cellIndex">
                <template v-if="typeof cell === 'string' && cell.startsWith('data:image')">
                  <img :src="cell" />
                </template>
                <template v-else>
                  {{ cell }}
                </template>
              </td>
            </tr>
          </tbody>
        </table>
      </div>
    </div>
  </template>
  
  <script>
  import axios from 'axios';
  
  export default {
    data() {
      return {
        activeTab: 'smiles',
        reactants: '',
        taskName: '',
        reactionType: '',
        results: null,     // 保存预测结果
        dataFrame: null    // 保存DataFrame结果
      };
    },
    methods: {
      async submitForm() {
        try {
          // 发送请求到后端
          const response = await axios.post('http://localhost:5002/predict', {
            reactants: this.reactants,
            taskName: this.taskName,
            reactionType: this.reactionType,
          });
  
          // 打印响应，确认数据结构
          console.log(response.data);
  
          // 解析返回的结果并显示在前端
          this.displayResults(response.data.results_dict);
          this.displayDataFrame(response.data.results_df);
        } catch (error) {
          console.error("Error submitting form:", error);
        }
      },
      displayResults(data) {
        // 将返回的结果存储到 `results` 中，Vue 会自动更新DOM
        this.results = data;
      },
      displayDataFrame(df) {
        // 将返回的DataFrame JSON字符串解析为对象
        this.dataFrame = JSON.parse(df);
      }
    }
  };
  </script>
  
  <style scoped>
  .container {
    width: 100%;
    margin: 0 auto;
    padding: 20px;
  }
  
  .header {
    background-color: #4086f4;
    color: white;
    padding: 10px;
    text-align: center;
  }
  
  .nav {
    display: flex;
    justify-content: space-between;
    background-color: #dcdcdc;
    padding: 10px;
  }
  
  .nav a {
    text-decoration: none;
    color: black;
    padding: 10px;
  }
  
  .content {
    display: flex;
    margin-top: 20px;
  }
  
  .form-container {
    flex: 1;
    padding: 20px;
    border: 1px solid #dcdcdc;
    margin-right: 20px;
  }
  
  .form-container label,
  .form-container input {
    display: block;
    width: 100%;
    margin-bottom: 10px;
  }
  
  .form-container input[type="text"] {
    padding: 10px;
    border: 1px solid #ccc;
    border-radius: 4px;
  }
  
  .form-container input[type="submit"] {
    background-color: #4086f4;
    color: white;
    border: none;
    padding: 10px;
    border-radius: 4px;
    cursor: pointer;
  }
  
  .tab-container {
    display: flex;
    justify-content: space-between;
    margin-bottom: 20px;
  }
  
  .tab-container button {
    flex: 1;
    padding: 10px;
    border: 1px solid #dcdcdc;
    background-color: #f0f0f0;
    cursor: pointer;
  }
  
  .tab-container button.active {
    background-color: #e0e0e0;
  }
  
  table {
    width: 100%;
    border-collapse: collapse;
  }
  
  table, th, td {
    border: 1px solid black;
  }
  
  th, td {
    padding: 10px;
    text-align: left;
  }
  
  img {
    max-width: 100px;
    height: auto;
  }
  </style>
  