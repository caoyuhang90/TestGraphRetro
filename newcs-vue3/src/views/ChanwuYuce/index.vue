<template>
  <div class="container">
    <div class="header">
      <h1>反应产物结构预测</h1>
    </div>

    <div class="nav">
      <router-link to="#">创建任务</router-link>
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
          <label for="reactant-count">输入反应物数量（小于50）</label>
          <input type="number" v-model="reactantCount" id="reactant-count" min="1" max="49" required />

          <div v-for="index in reactantCount" :key="'reactant-' + index">
            <label :for="'reactant-' + index">反应物 {{ index }}</label>
            <input type="text" v-model="reactants[index - 1]" :id="'reactant-' + index" required />
          </div>

          <label for="solvent-count">输入溶剂数量（小于50）</label>
          <input type="number" v-model="solventCount" id="solvent-count" min="1" max="49" required />

          <div v-for="index in solventCount" :key="'solvent-' + index">
            <label :for="'solvent-' + index">溶剂 {{ index }}</label>
            <input type="text" v-model="solvents[index - 1]" :id="'solvent-' + index" required />
          </div>

          <label for="task-name">任务名称</label>
          <input type="text" v-model="taskName" id="task-name" />

          <label for="reaction-type">反应类型（选填）</label>
          <input type="text" v-model="reactionType" id="reaction-type" />

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
      reactantCount: 1,
      solventCount: 1,
      reactants: Array(49).fill(''), // 预定义最大数量的反应物
      solvents: Array(49).fill(''),   // 预定义最大数量的溶剂
      taskName: '',
      reactionType: '',
      results: null,
      dataFrame: null,
    };
  },
  methods: {
    async submitForm() {
      try {
        const combinedSmiles = [...this.reactants, ...this.solvents].filter(Boolean).join('.');
        
        const response = await axios.post('http://localhost:5002/predict', {
          reactants: combinedSmiles,
          taskName: this.taskName,
          reactionType: this.reactionType,
        });

        console.log(response.data);
        this.displayResults(response.data.results_dict);
        this.displayDataFrame(response.data.results_df);
      } catch (error) {
        console.error("Error submitting form:", error);
      }
    },
    displayResults(data) {
      this.results = data;
    },
    displayDataFrame(df) {
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
    background-color:rgb(1,70,138);
    color: white;
    padding: 10px;
    text-align: center;
  }
  
  .nav {
    display: flex;
    justify-content: space-between;
    background-color: #fff;
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
    background-color: #fff;
  }
  
  .result-container{
    background-color: #fff;
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
    width: 5%;
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
    background-color: #4086f4;
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
